# This pipeline simulates the out of africa model and finds AIMs using naiveBayes
masterdir = "./"
source(paste(masterdir, "code/tools.r", sep=""))
require("parallel")

setcores(cores = 3)

seed = 11 # in order to have reproducable results
testindsperdeme = 200
demes=3
informativenessBound = 0.9
stepsize = 10000
numAIM = 15 # number of AIMs 
method = "naiveBayes"
error.type="logloss"
diploid = TRUE
set.seed(seed)
size = rep(testindsperdeme, demes) # size of the test set in all demes

simulate = TRUE # Simulate datasets
prepare = TRUE # Make Training and Test sets
step1 = TRUE # Get frequencies for all SNPs in all populations for the training set
step2 = TRUE # Filter most informative SNPs and write a new vcf File
step3 = TRUE # Cut files in smaller pieces for better handling
step4 = TRUE # Get nss for all chromosomes
step5 = TRUE # Store training data in RDS-Files
step6 = TRUE  # Read data in chunks and find possible AIMs
step7 = TRUE  # Create a vcf.gz File for all SNPs from step6; this is called AIMs.vcf.gz
step8 = TRUE # Create classification and posterior probabilities


if(simulate){
  print("simulating SNP data for out of africa model -- this may take a while...")
  globalcommand = "python3 code/simulate4biogeo.py ooa"
  commands = c()
  allseeds = seed
  # to simulate all seeds used for the publication at once uncomment:
  # allseeds = c(11,12,13,14,15,16,17,18,19,20)
  for(itseed in allseeds){
    commands = c(commands,paste(globalcommand, itseed, 20,  sep= " "))
  }
  exec<-function(commandline){
    system(commandline)
  }
  library(parallel)
  mclapply(commands,exec)
  for(itseed in seeds){
    command = paste("mv ooa_chromosome_*_seed_",itseed,".vcf.gz data/sim/ooa/tmp/.",sep="")
    system(command)
  }
}


# Make Training and Test sets
filenameInds = "ooa_SampleListWithLocations.txt"
filenameIndsTraining = "training_ooa_SampleListWithLocations.txt"
filenameIndsTest = "test_ooa_SampleListWithLocations.txt"

filenameData = NULL
for(i in 1:20) {
  filenameData = c(filenameData, paste("data/sim/ooa/tmp/ooa_chromosome_", i, "_seed_", seed, sep=""))
}

if(prepare) {
  a = read.vcfR(paste("data/sim/ooa/tmp/ooa_chromosome_1_seed_", seed, ".vcf.gz", sep=""), nrows=1)
  inds = colnames(a@gt)[-1]
  inds = cbind(inds, rep(0, 2400), c(rep("AFR", 800), rep("EUR", 800), rep("ASI", 800)))
  write.table(inds, filenameInds, col.names= FALSE, row.names = FALSE)
  getTestandtrainingset(filenameInds, filenameIndsTraining, filenameIndsTest, size) 
}

# load frequencies in all pops into R

if(step1) {   # Get frequencies in all superpopulations
  for(i in 20:1) {
    getfreqs.from.vcf(filenameData[i], filenameIndsTraining, snplist=FALSE) 
  }
  # mclapply(filenameData,FUN = function(x){getfreqs.from.vcf(x,filenameIndsTraining,snplist = FALSE)})  # this does not work yet
}

# get positions of most informative SNPs and create new vcf.gz files

if(step2) {
  #  source("../biogeoanc-prediction/util.r")
  for(i in 20:1) {
    load(paste(filenameData[i], "_freqs", sep=""))
    info = informativeness.global(freqs)
    nss = ncol(freqs)
    
    informativenessQuantile = sort(info, decreasing = TRUE)[(1-informativenessBound) * nss]
    chr.informative = chr[info>informativenessQuantile]
    pos.informative = pos[info>informativenessQuantile]
    options(scipen=999) # This prevents 48000000 to be written as 4.8e+07 which cannot be read again
    write.table(cbind(chr.informative, pos.informative), file = paste("tmp/chr", i, "_informativePositions", sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE)
    options(scipen=0) 
    
    call = paste("vcftools --gzvcf ",filenameData[i],".vcf.gz --chr ", i-1, " --positions tmp/chr", i, "_informativePositions --recode-INFO-all --recode --out tmp/chr", i, "_seed_", seed, "_informative", sep="")
    cat("Getting most informative positions on chromosome ", i, "\n")
    system(call)
    cat("Zipping file for chromosome ", i, "\n")
    call=paste("gzip -f tmp/chr", i, "_seed_", seed, "_informative.recode.vcf", sep="")
    system(call)
  }
}

# Cut files in smaller pieces and get rid of unnamed SNPs
# This results in tmp/chr1_segment0.recode.vcf.gz etc...

if(step3) {
  filenames = NULL
  for(i in 1:20) {
    filenames = c(filenames, paste("tmp/chr", i, "_seed_", seed, "_informative.recode.vcf.gz", sep=""))
  }
  cut.VCF(filenames, stepsize, removeUnnamed = TRUE)
}

# get file nss, which stores the number of informative segregating sites for all chromosomes

if(step4) {
  nss = NULL
  for(i in 1:20) {
    nss = c(nss, getNss(paste("tmp/chr", i, "_seed_", seed, "_informative.recode.vcf.gz", sep="")))
  }
  write.table(nss, file = paste("nss_seed_", seed, sep=""), quote=FALSE, col.names=FALSE, row.names=FALSE)
}

# store data for trainingset in RDS-Files

if(step5) {
  a = read.csv(filenameIndsTraining, sep="", header=FALSE)
  training.inds = a[,c(1,3)]
  nss = read.table(paste("nss_seed_", seed, sep=""))[,1]
  
  filenames = NULL
  for(i in 1:20) {
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
  a = read.csv(filenameIndsTraining, sep="", header=FALSE)
  training.inds = a[,c(1,3)]
  nss = read.table(paste("nss_seed_", seed, sep=""))[,1]
  
  filenames = NULL
  for(i in 1:20) {
    for(j in 0:((nss[i]-1)/stepsize)) {
      filenames = c(filenames, paste("tmp/chr", i, "_segment", j, ".recode.vcf.gz.rds", sep=""))
    }
  }
  sampleAIMs = NULL
  sample.AIMs = find.AIMs.fromRDS(filenames, sampleAIMs, training.inds, numAIM, method, error.type, diploid)
}

# We create a vcf.gz File for all SNPs from step6; this is called AIMs.vcf.gz
if(step7) {
  filenameDataVCFGZ = filenameData
  for(i in 1:20) filenameDataVCFGZ[i] = paste(filenameData[i], ".vcf.gz", sep="")
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
