masterdir = "./"
source(paste(masterdir, "code/tools.r", sep=""))
require("parallel")

setcores(cores = 3)

seed = 11 # in order to have reproducable results
testindsperdeme = 20
demes=3
informativenessBound = 0.9
stepsize = 10000
numAIM = 10 # number of AIMs 
method = "naiveBayes"
error.type="logloss"
diploid = TRUE



set.seed(seed)
size = rep(testindsperdeme, demes) # size of the test set in all demes

prepare0 = TRUE # Make Training and Test sets
step1 = TRUE # Get frequencies for all SNPs in all populations for the training set
step2 = TRUE # Filter most informative SNPs and write a new vcf File
step3 = TRUE # Cut files in smaller pieces for better handling
step4 = TRUE # Get nss for all chromosomes
step5 = TRUE # Store training data in RDS-Files
step6 = TRUE  # Read data in chunks and find possible AIMs
step7 = TRUE  # Create a vcf.gz File for all SNPs from step6; this is called AIMs.vcf.gz
step8 = TRUE # Create classification and posterior probabilities



# Make Training and Test sets
filenameInds = "ooa_SampleListWithLocations.txt"
filenameIndsTraining = "training_ooa_SampleListWithLocations.txt"
filenameIndsTest = "test_ooa_SampleListWithLocations.txt"

filenameData = "data/sim/ooa/ooa_chromosome_1_example"

if(prepare0) {
  a = read.vcfR(paste(filenameData,".vcf.gz",sep=""), nrows=1)
  inds = colnames(a@gt)[-1]
  inds = cbind(inds, rep(0, 240), c(rep("AFR", 80), rep("EUR", 80), rep("ASI", 80)))
  write.table(inds, filenameInds, col.names= FALSE, row.names = FALSE)
  getTestandtrainingset(filenameInds, filenameIndsTraining, filenameIndsTest, size) 
}

# load frequencies in all pops into R

if(step1) {   # Get frequencies in all superpopulations
    getfreqs.from.vcf(filenameData, filenameIndsTraining, snplist=FALSE) 
}

# get positions of most informative SNPs and create new vcf.gz files

if(step2) {
    load(paste(filenameData, "_freqs", sep=""))
    info = informativeness.global(freqs)
    nss = ncol(freqs)

    informativenessQuantile = sort(info, decreasing = TRUE)[(1-informativenessBound) * nss]
    chr.informative = chr[info>informativenessQuantile]
    pos.informative = pos[info>informativenessQuantile]
    options(scipen=999) # This prevents 48000000 to be written as 4.8e+07 which cannot be read again
    write.table(cbind(chr.informative, pos.informative), file = "tmp/chr1_informativePositions", quote=FALSE, col.names=FALSE, row.names=FALSE)
    options(scipen=0) 

    call = paste("vcftools --gzvcf ", filenameData, ".vcf.gz --positions tmp/chr1_informativePositions --recode-INFO-all --recode --out tmp/chr1_informative", sep="")
    cat("Getting most informative positions on chromosome ", 1, "\n")
    system(call)
    cat("Zipping file for chromosome ", 1, "\n")
    call=paste("gzip -f tmp/chr1_informative.recode.vcf", sep="")
    system(call)
}

# Cut files in smaller pieces and get rid of unnamed SNPs
# This results in tmp/chr1_segment0.recode.vcf.gz etc...

if(step3) {
  filenames = paste("tmp/chr1_informative.recode.vcf.gz", sep="")
  cut.VCF(filenames, stepsize, removeUnnamed = TRUE)
}

# get file nss, which stores the number of informative segregating sites for all chromosomes

if(step4) {
  nss = getNss("tmp/chr1_informative.recode.vcf.gz")
  write.table(nss, file = "nss.tab", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# store data for trainingset in RDS-Files

if(step5) {
  a = read.csv(filenameIndsTraining, sep="", header=FALSE)
  training.inds = a[,c(1,3)]
  nss = read.table("nss.tab")[,1]
    filenames = NULL
    for(j in 0:((nss-1)/stepsize)) {
      filenames = c(filenames, paste("tmp/chr1_segment", j, ".recode.vcf.gz", sep=""))
    }
  for(i in 1:length(filenames)) {
    dat = getData(filenames[i], inds = training.inds, diploid = TRUE)
    saveRDS(dat$sample, file = paste(filenames[i], ".rds", sep=""), version = 2)
  }
}

# dat.training contains data for the training set. The order of rows is the same as in filenameIndsTraining

if(step6) {
  a = read.csv(filenameIndsTraining, sep="", header=FALSE)
  training.inds = a[,c(1,3)]
  nss = read.table("nss.tab")[,1]

  filenames = NULL
  for(j in 0:((nss[i]-1)/stepsize)) {
      filenames = c(filenames, paste("tmp/chr1_segment", j, ".recode.vcf.gz.rds", sep=""))
    }
  sampleAIMs = NULL
  sample.AIMs = find.AIMs.fromRDS(filenames, sampleAIMs, training.inds, numAIM, method, error.type, diploid)
}

# We create a vcf.gz File for all SNPs from step6; this is called AIMs.vcf.gz
if(step7) {
  filenameDataVCFGZ = paste(filenameData, ".vcf.gz", sep="")
  extract.and.merge(filenameDataVCFGZ, "AIMs", outfile="AIMs.vcf")
  system("gzip -f AIMs.vcf")
}

# compute posterior probabilities with final AIMset and classify samples with naive Bayes classifier
if(step8) {
  a = read.csv(filenameInds, sep=" ", header=FALSE)
  inds = a[,c(1,2,3)]
  samLoc = a[,3]
  a = getData("AIMs.vcf.gz", skip=0, windowSize=-1, inds = inds, snps = TRUE, diploid = TRUE)
  sample = a$sample
  freqs = getfreqs(sample, samLoc, diploid = TRUE)
  samplesInDeme = table(samLoc)
  predictions = prediction.naiveBayes(sample, freqs, samplesInDeme, diploid = TRUE)
  rownames(predictions) = inds[,1]
  colnames(predictions) = names(table(samLoc))
  classifications = rownames(freqs)[apply(predictions, 1, which.max)]
  names(classifications) = inds[,1]
  write.csv(predictions, file = "predictions.csv" )
  write.table(classifications, file = "classifications.tab", col.names = F )
}

