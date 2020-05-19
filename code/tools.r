require("vcfR")
require("parallel")

# Approximation of the test error

appTestError<-function(sample.test, freqs.training, samplesInDeme.training, samLoc.test, method = "naiveBayes", error.type="misclassification", diploid = FALSE){

  if(method == "naiveBayes") {
    predictions = prediction.naiveBayes(sample.test, freqs.training, samplesInDeme.training, diploid)
    colnames(predictions) = sort(unique(samLoc.test))
  }
  if(error.type=="misclassification") {
    error = misclassification.error(predictions, samLoc.test)
  }
  if(error.type=="logloss") {
#    error = logloss.error(predictions, as.numeric(samLoc.test))
    error = logloss.error(predictions, samLoc.test)
  }
  error
}

# It is often easier to handle datafiles when they are cut into smaller files;

cut.VCF<-function(filenamesDataVCFGZ, stepsize, removeUnnamed = TRUE) {
  for(i in 1:length(filenamesDataVCFGZ)) {
    cat("Cutting file ", filenamesDataVCFGZ[i], ".\n")
    identifiers = getIdentifiers(filenamesDataVCFGZ[i])
    if(removeUnnamed) {
      # get rid of unnamed SNPs
      identifiers = identifiers[identifiers != "."]
    }
    
    nss = length(identifiers)
    for(j in 0:((nss-1)/stepsize)) {
      options(scipen=999) # This prevents 48000000 to be written as 4.8e+07 which cannot be read again
      write.table(identifiers[j * stepsize + 1:(min(stepsize, nss - j*stepsize))], file = "next.identifiers", quote=FALSE, col.names=FALSE, row.names=FALSE)
      options(scipen=0)
      cat("Writing segment ", j, "(from", as.integer((nss-1)/stepsize), ").\n")

      call = paste("vcftools --gzvcf ", filenamesDataVCFGZ[i], " --snps next.identifiers --recode --out ", "tmp/chr", i, "_segment", j, sep="")
      system(call)
      system(paste("gzip -f tmp/chr", i, "_segment", j, ".recode.vcf", sep=""))   
    }
  }
  system("rm next.identifiers")
}

# From all files of the vector filenames we extract the SNPs in snplist from the individuals in inds,
# and merge to a single file
# all this is done using vcftools and bcftools

extract.and.merge<-function(filenamesDataVCFGZ, snplist, outfile) {
  for(i in 1:length(filenamesDataVCFGZ)) {
    call = paste("vcftools --gzvcf ", filenamesDataVCFGZ[i], " --snps ", snplist, " --recode --out ", filenamesDataVCFGZ[i], "_snps_extracted", sep="")
    system(call)
  }
  call = paste("bcftools concat -Ov ")
  for(i in 1:length(filenamesDataVCFGZ)) {
    call = paste(call, paste(filenamesDataVCFGZ[i], "_snps_extracted.recode.vcf", sep=""))
  }
  call = paste(call, " -o ", outfile)
  system(call)
  for(i in 1:length(filenamesDataVCFGZ)) {
    system(paste("rm -r ", filenamesDataVCFGZ[i], "_snps_extracted.recode.vcf", sep=""))
  }
}

# filenames is a vector of RDS-files, which contain all samples for the training individuals

find.AIMs.fromRDS<-function(filenamesRDS, sampleAIMs = NULL, inds, numAIM, method = "naiveBayes", error.type="misclassification", diploid = FALSE, outfile = "AIMs"){

  samLoc = inds[,2]
  samplesInDeme = table(samLoc)
  
  for(i in 1:numAIM) {
    best = mclapply(1:length(filenamesRDS),function(x){
        find.nextAIM.fromRDS(filenamesRDS[x], sampleAIMs, inds, method, error.type, diploid)})
    error = rep(0, length(filenamesRDS))
    for(j in 1:length(error)) error[j] = best[[j]]$min.error
    bestSNP = best[[which.min(error)]]$bestSNP
    min.error = min(error)
    sample.nextAIM = best[[which.min(error)]]$sample.nextAIM
    sampleAIMs = cbind(sampleAIMs, sample.nextAIM)
    AIMs = colnames(sampleAIMs)
    options(scipen=999) # This prevents 48000000 to be written as 4.8e+07 which cannot be read again
    write.table(AIMs, file = outfile, quote=FALSE, col.names=FALSE, row.names=FALSE)
    options(scipen=0)
    saveRDS(sampleAIMs, file = "sampleAIMs.rds", version = 2)    
  }
  sampleAIMs
}

# inds contains the training set, filename is a RDS-file

find.nextAIM.fromRDS<-function(filenameRDS, sampleAIMs, inds, method, error.type, diploid) {
  samLoc = inds[,2]
  samplesInDeme = table(samLoc)
  no = if(is.null(ncol(sampleAIMs))) 0 else ncol(sampleAIMs)
  min.error = 1000
  
  cat("Filename: ", filenameRDS, "\n")
  cat("We already have ", no, " AIMs.\n")
  sample = readRDS(filenameRDS)
  sampleAIMandSNPs = cbind(sampleAIMs, sample)
  freqs.training = getfreqs(sampleAIMandSNPs, samLoc, diploid)

  if(is.null(sampleAIMs)) {
    error<-unlist(lapply(1:ncol(sample),function(x){
        appTestError(sampleAIMandSNPs[,x], freqs.training[,x], samplesInDeme, samLoc, method, error.type, diploid)}))
  } else {
    error<-unlist(lapply(no + (1:ncol(sample)),function(x){
        appTestError(sampleAIMandSNPs[,c(1:no, x)], freqs.training[,c(1:no, x)], samplesInDeme, samLoc, method, error.type, diploid)}))
  }
#  gc()
#  names(error) = colnames(sample)

  min.error = min(error)
  bestSNP = which.min(error)
  sample.nextAIM = as.matrix(sample[,bestSNP], ncol=1)
  colnames(sample.nextAIM) = colnames(sample)[bestSNP]
  rownames(sample.nextAIM) = inds[,1]
  
  cat("In file ", filenameRDS, " the best AIM is ", colnames(sample)[bestSNP], ". The training error is ", min.error, ".\n", sep="")

  cat("The misclassification error on the training set is ", appTestError(cbind(sampleAIMs, sample.nextAIM), getfreqs(cbind(sampleAIMs, sample.nextAIM), samLoc, diploid), samplesInDeme, samLoc, method, error.type = "misclassification", diploid), ".\n\n")

  list(filename = filenameRDS, bestSNP = colnames(sample)[bestSNP], min.error = min.error, sample.nextAIM = sample.nextAIM)
}

# Determine which SNPs from a snplist occur in a set of files filenamesDataVCFGZ;

find.SNPs.inVCFGZ<-function(filenamesDataVCFGZ, snplist) {
  SNPs = as.vector(read.table(snplist, header=FALSE)[,1])
  jointSNPs = NULL
  inSNPlistNotInData = NULL
  for(i in 1:length(filenamesDataVCFGZ)) {
    cat("Finding joint SNPs in file ", filenamesDataVCFGZ[i], "\n")
    identifiers = as.vector(getIdentifiers(filenamesDataVCFGZ[i]) )
    jointSNPs = c(jointSNPs, intersect(identifiers, SNPs))
  }
  list(jointSNPs = jointSNPs, inSNPlistNotInData = setdiff(SNPs, jointSNPs))  
}

# inds is nsam x 2 matrix, where the first column are the column names from the vcf.gz file which we want to read;
# the second column is the sampling location
# snps is a vector of identifiers which SNPs should come back;

getData<-function(filenameDataVCFGZ, skip=0, windowSize=-1, inds, snps = TRUE, diploid = FALSE) {
  if(diploid) {
    trans = c(0,1,1,2, NA)
    names(trans)=c("0|0", "0|1", "1|0", "1|1", "NA")
  } else {
    trans = c(0,1, NA)
    names(trans)=c("0", "1", "NA")
  }
  samLoc = inds[,2]
  
  a = read.vcfR(filenameDataVCFGZ, skip = skip, nrows = windowSize)
  genotypes = a@gt[,-1]
  row.names(genotypes) = a@fix[,3]
  sample = matrix(0, ncol = nrow(genotypes), nrow = nrow(inds))
  rownames(sample) = as.character(inds[,1])
  colnames(sample) = as.character(a@fix[,3])
  for(k in 1:nrow(inds)) {
    sample[k,] = as.numeric(trans[genotypes[,as.character(inds[k,1])]])
  }
  if(snps) {
    snps = colnames(sample) 
    for(i in c(which(snps == ""),which(is.na(snps)))) {
      colnames(sample)[i] = snps[i] = paste(filenameDataVCFGZ, "_", skip + i, sep="") # sometimes SNPs have no identifier
    }
  }  
  res = list(sample = sample[, snps], samLoc = samLoc)
  res
}

# function to get the frequencies of SNPs per deme
getfreqs = function(sample, samLoc, diploid = FALSE){
  if(ncol(as.matrix(sample))==1) {
    as.matrix(by(as.matrix(sample),samLoc,colMeans)) / (1+diploid)
  } else {
    t(simplify2array(by(sample,samLoc,colMeans))) / (1+diploid)
  }  
}

# same as last function, but vcftoolsis used to compute frequencies. This function is
# used on large datafiles.

getfreqs.from.vcf<-function(filenameData, filenameInds, snplist=FALSE) {
  a = as.matrix(read.table(filenameInds, header=FALSE))
  pops = sort(unique(a[,3]))
  inds = a[,c(1,3)]
  noinds = rep(0, length(pops)) # will contain the number of samples per population

  for(i in 1:length(pops)) {
    call = paste("vcftools --gzvcf ", filenameData, ".vcf.gz", sep="")
    if(snplist) call = paste(call, " --snps ", snplist)
    indslocal = inds[inds[,2]==pops[i],1]  
    noinds[i] = length(indslocal)
    for(j in 1:length(indslocal)) {
      call = paste(call, " --indv ", indslocal[j])
    }
    call = paste(call, " --freq2 --out ", filenameData, "_", pops[i], sep="")
    system(call)
  }

  a<-read.csv(paste(filenameData, "_", pops[1], ".frq", sep=""), header=TRUE, sep="\t", row.names=NULL)
  chr = as.numeric(a[,1])
  pos = as.numeric(a[,2])
  freqs = matrix(0, nrow = length(pops), ncol = length(chr))
  for(i in 1:length(pops)) {
    cat("Storing ", paste(filenameData, "_", pops[i], ".frq", sep=""), " in R format\n")
    a<-read.csv(paste(filenameData, "_", pops[i], ".frq", sep=""), header=TRUE, sep="\t", row.names=NULL)
    freqs[i,] = as.numeric(a[,5])
  }
  cat("Storing chr, pos, pops, freqs, noinds in file ", filenameData, "_freqs\n")
  save(chr, pos, pops, freqs, noinds, file = paste(filenameData, "_freqs", sep=""))
}

# Obtain the list of SNP identifiers from a file;

getIdentifiers<-function(filenameDataVCFGZ) {
  call = paste("gunzip -c ", filenameDataVCFGZ, " | grep \"^[^#]\" | cut -f 3 > identifiers.out")
  system(call)
  identifiers = read.table("identifiers.out", header=FALSE)
  system("rm identifiers.out")
  unlist(identifiers)
}

# Obtains the number of segregating sites from a file; 
getNss<-function(filenameDataVCFGZ) {
  call = paste("gunzip -c ", filenameDataVCFGZ, " | grep \"^[^#]\" | wc -l > nss.out")
  system(call)
  nss = read.table("nss.out")
  system("rm nss.out")
  unlist(nss)
}

# This function outputs two new files called filenameIndsTraining and filenameIndsTest.
# size is a vector containing the sizes of the testset in all populations

getTestandtrainingset<-function(filenameInds, filenameIndsTraining, filenameIndsTest, size) {
  a = as.matrix(read.table(filenameInds, header=FALSE))
  pops = sort(unique(a[,3]))
  inds = a[,c(1,3)]
  noinds = rep(0, length(pops)) # will contain the number of samples per population
  testset=trainingset = NULL
  for(i in 1:length(pops)) {
    indslocal = inds[inds[,2]==pops[i],1]  
    noinds[i] = sum(inds[,2]==pops[i])
    s = sample(1:noinds[i], size[i], replace = FALSE)
    testset = c(testset, indslocal[s])
    trainingset = c(trainingset, indslocal[-s])
  }

  system(paste("rm ", filenameIndsTest))
  system(paste("rm ", filenameIndsTraining))
  for(i in 1:length(testset)) {
    call = paste("grep -w ", testset[i], filenameInds, " >> ", filenameIndsTest)
    system(call)
  }
  for(i in 1:length(trainingset)) {
    call = paste("grep -w ", trainingset[i], filenameInds, " >> ", filenameIndsTraining)
    system(call)
  }
}

# This function can handle the result of getfreqs 
# it gives the informativeness for all SNPs
informativeness.global = function(freqs){
  mfreqs = colMeans(freqs)
  as.vector(- xlogx2(mfreqs) + colMeans(xlogx2(freqs)))
}

# here is the logloss error 
logloss.error<-function(prediction, true.value) {
  # these are the probabilites for the correct classes
  probs.for.correct.classes<-diag(prediction[,true.value])
  -mean(log(probs.for.correct.classes))
}

# Here is the light version of STRUCTURE, similar to Snipper
# sample.test is a matrix
# We use prediction.naiveBayes(...) = exp(log.prediction.naiveBayes(...)) for numerical stability

log.prediction.naiveBayes<-function(sample.test, freqs.training, samplesInDeme.training, diploid) {
  demes = nrow(as.matrix(freqs.training))
  counts = freqs.training * as.vector(samplesInDeme.training) * (1+diploid) # counts now contains natural numbers
  
  res = matrix(0, nrow = nrow(as.matrix(sample.test)), ncol = demes)

  for(k in 1:demes) {
    # This is an implementation of formula (1)	 in the main text
    res[,k] = apply(as.matrix(log((1 + t(counts)[,k])) * t(sample.test) + log((1+diploid) * samplesInDeme.training[k] - t(counts)[,k]+1)*(1 + diploid -t(sample.test))), 2, sum)
    # res[,k] = apply(as.matrix(log((1 + t(counts)[,k]) * t(sample.test) +
    #	    ((1+diploid) * samplesInDeme.training[k] - t(counts)[,k]+1)*(1 + diploid -t(sample.test)))), 2, sum)
  }
  res 
}

# here is the misclassification error 
misclassification.error<-function(prediction, true.value) {
  mean(colnames(prediction)[apply(prediction, 1, which.max)] != true.value)
}

# See log.prediction;

prediction.naiveBayes<-function(sample.test, freqs.training, samplesInDeme.training, diploid) {
  res = log.prediction.naiveBayes(sample.test, freqs.training, samplesInDeme.training, diploid)
  res = exp(res - max(res))
  res / apply(res, 1, sum)
}

# Since finding AIMs has multicore support, it is advisable to set the number of cores to a maximum.

setcores <- function(cores=1) {
  options(mc.cores = cores) # how many cores should be used 
  if(getOption("mc.cores") < detectCores()-1)
    warning("This system provides more cores than used, consider using more cores to speed up the search")
  if(getOption("mc.cores") == detectCores())
    warning("You are using all available cores. Maybe consider using one core less.")
  if(getOption("mc.cores") > detectCores())
    stop("You are using more cores than available on this machine. Lower the number of used cores.")
}

# The following two functions are used in informativeness
xlogx2<-function(x) {
  xlogx(x) + xlogx(1-x)
} 
xlogx<-function(x) {
  ifelse(x==0, 0, x *log(x))
} 

