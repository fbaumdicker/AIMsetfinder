masterdir = "./"
source(paste(masterdir, "code/tools.r", sep=""))
require("parallel")

setcores(cores = 3)

set.seed(31) # in order to have reproducible results
testindsperdeme = 100
demes=4
informativenessBound = 0.9
stepsize = 10000
numAIM = 30 # number of AIMs 
method = "naiveBayes"
error.type="logloss"
diploid = TRUE

size = rep(testindsperdeme, demes) # size of the test set in all demes

prepare = TRUE # Get rid of AMR-individuals, make Training and Test sets
step0 = TRUE # reduce dataset to bi-allelic SNPs, and remove indels
step1 = TRUE # Get frequencies for all SNPs in all populations for the training set
step2 = TRUE # Filter most informative SNPs and write a new vcf File
step3 = TRUE # Cut files in smaller pieces for better handling
step4 = TRUE # Get nss for all chromosomes
step5 = TRUE # Store training data in RDS-Files
step6 = TRUE  # Read data in chunks and find possible AIMs
step7 = TRUE  # Create a vcf.gz File for all SNPs from step6; this is called AIMs.vcf.gz
step8 = TRUE # Create classification and posterior probabilities

filenameInds = "1000G_SampleListWithLocations_noAMR.txt"
filenameIndsTraining = "training_1000G_SampleListWithLocations_noAMR.txt"
filenameIndsTest = "test_1000G_SampleListWithLocations_noAMR.txt"
filenameData = NULL
for(i in 1:22) {
  filenameData = c(filenameData, paste("tmp/chr", i, "_biallelic.recode", sep=""))
}

# Get rid of AMR-individuals, make Training and Test sets
if(prepare) {
  system("grep -E 'AFR|EUR|SAS|EAS' data/1000G/1000G_SampleListWithLocations.txt > 1000G_SampleListWithLocations_noAMR.txt")
  # This stores the training and test sets
  getTestandtrainingset(filenameInds, filenameIndsTraining, filenameIndsTest, size) 
}

# This results in files tmp/chr22_biallelic.vcf.gz etc.
if(step0) {
  for(i in 22:1) {
    call = paste("vcftools --gzvcf ", masterdir, "data/1000G/ALL.chr", i, ".*.vcf.gz --min-alleles 2 --max-alleles 2 --remove-indels --recode-INFO-all --recode --out tmp/chr", i, "_biallelic", sep="")
    system(call)
    call=paste("gzip -f tmp/chr", i, "_biallelic.recode.vcf", sep="")
    system(call)
  }
}

# load frequencies in all pops into R

if(step1) {   # Get frequencies in all superpopulations
  for(i in 22:1) {
    getfreqs.from.vcf(filenameData[i], filenameIndsTraining, snplist=FALSE) 
  }
}

# get positions of most informative SNPs and create new vcf.gz files

if(step2) {
#  source("../biogeoanc-prediction/util.r")
  for(i in 22:1) {
    load(paste(filenameData[i], "_freqs", sep=""))
    info = informativeness.global(freqs)
    nss = ncol(freqs)

    informativenessQuantile = sort(info, decreasing = TRUE)[(1-informativenessBound) * nss]
    chr.informative = chr[info>informativenessQuantile]
    pos.informative = pos[info>informativenessQuantile]
    options(scipen=999) # This prevents 48000000 to be written as 4.8e+07 which cannot be read again
    write.table(cbind(chr.informative, pos.informative), file = paste("tmp/chr", i, "_informativePositions", sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE)
    options(scipen=0) 

    call = paste("vcftools --gzvcf tmp/chr", i, "_biallelic.recode.vcf.gz --chr ", i, " --positions tmp/chr", i, "_informativePositions --recode-INFO-all --recode --out tmp/chr", i, "_informative", sep="")
    cat("Getting most informative positions on chromosome ", i, "\n")
    system(call)
    cat("Zipping file for chromosome ", i, "\n")
    call=paste("gzip -f tmp/chr", i, "_informative.recode.vcf", sep="")
    system(call)
  }
}

# Cut files in smaller pieces and get rid of unnamed SNPs
# This results in tmp/chr1_segment0.recode.vcf.gz etc...

if(step3) {
  filenames = NULL
  for(i in 1:22) {
    filenames = c(filenames, paste("tmp/chr", i, "_informative.recode.vcf.gz", sep=""))
  }
  cut.VCF(filenames, stepsize, removeUnnamed = TRUE)
}

# get file nss, which stores the number of informative segregating files for all chromosomes

if(step4) {
  nss = NULL
  for(i in 1:22) {
    nss = c(nss, getNss(paste("tmp/chr", i, "_informative.recode.vcf.gz", sep="")))
  }
  write.table(nss, file = "nss", quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# store data for trainingset in RDS-Files

if(step5) {
  a = read.csv(filenameIndsTraining, sep="\t", header=FALSE)
  training.inds = a[,c(1,3)]
  nss = read.table("nss")[,1]

  filenames = NULL
  for(i in 1:22) {
    for(j in 0:((nss[i]-1)/stepsize)) {
      filenames = c(filenames, paste("tmp/chr", i, "_segment", j, ".recode.vcf.gz", sep=""))
    }
  }
  for(i in 1:length(filenames)) {
    dat = getData(filenames[i], inds = training.inds, diploid = TRUE)
    saveRDS(dat$sample, file = paste(filenames[i], ".rds", sep=""), version = 2)
  }
}

# dat.training contains data for the training set. The order of rows is the same as in filenameIndsTraining

if(step6) {
  a = read.csv(filenameIndsTraining, sep="\t", header=FALSE)
  training.inds = a[,c(1,3)]
  nss = read.table("nss")[,1]

  filenames = NULL
  for(i in 1:22) {
    for(j in 0:((nss[i]-1)/stepsize)) {
      filenames = c(filenames, paste("tmp/chr", i, "_segment", j, ".recode.vcf.gz.rds", sep=""))
    }
  }
  sampleAIMs = readRDS(file = "sampleAIMs.rds")
  #sampleAIMs = NULL
  sample.AIMs = find.AIMs.fromRDS(filenames, sampleAIMs, training.inds, numAIM, method, error.type, diploid)
}

# We create a vcf.gz File for all SNPs from step3; this is called possAIMs.vcf.gz
if(step7) {
  filenameDataVCFGZ = filenameData
  for(i in 1:22) filenameDataVCFGZ[i] = paste(filenameData[i], ".vcf.gz", sep="")
  extract.and.merge(filenameDataVCFGZ, "AIMs", filenameInds, outfile="AIMs.vcf")
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
  cat("pipeline has finished. AIM set, predictions, and classifications have been generated.\n")
}
