require("vcfR")
source("../../code/tools.r")
# read reference data from Snipper input file

a = read.csv("illumina_55.csv", sep="\t", stringsAsFactors = FALSE, header = FALSE)
training.samLoc = a[6:nrow(a),2]
names(training.samLoc) = a[6:nrow(a),3]
training.samLocFine = a[6:nrow(a),2]
names(training.samLocFine) = a[6:nrow(a),3]
training.sample = a[6:nrow(a), 4:ncol(a)]
rownames(training.sample) = a[6:nrow(a), 3]
colnames(training.sample) = a[1,4:ncol(a)]
training.sample[training.sample=="NN"] = NA

# Abbreviate superpops
trans = c("AMR", "AFR", "EUR", "EAS", "SAS", "MEA", "OCE", "SIB")
#names(trans) = c("America", "Africa", "Europe", "East Asia", "South Asia", "Middle East", "Oceania", "Siberia")
names(trans) = c("America", "Africa", "Europe", "East Asian", "South Asian", "Middle East", "Oceanian", "Siberia")
training.samLoc = trans[training.samLoc]
names(training.samLoc) = rownames(training.sample)

#read test data
a = read.csv("SNPstwoVPs.csv", stringsAsFactors = FALSE, header = FALSE)
test.sample = a[2:nrow(a), 2:ncol(a)]
rownames(test.sample) = a[2:nrow(a), 1]
colnames(test.sample) = a[1,2:ncol(a)]
test.sample[test.sample=="NN"] = NA

# check if all SNPs in the test.sample are available in the training.sample
for(i in colnames(test.sample)) {
  if(!is.element(i, colnames(training.sample))) cat("Warning: ", i, " not in the training data.\n", sep="")
}
training.sample = training.sample[,colnames(test.sample)]

# check if all SNPs are bi-allelic
# also define reference alleles
ref.alleles = 1:ncol(test.sample)
names(ref.alleles) = colnames(test.sample)
used.AIMs = colnames(test.sample)

for(AIM in colnames(test.sample)) {
  dat = c(as.vector(training.sample[,AIM]), as.vector(test.sample[,AIM]))
  dat = dat[!is.na(dat)]
  dat = unlist(strsplit(dat, ","))
  dat = unlist(strsplit(dat, ""))
  ref.alleles[AIM] = dat[1]
  if(length(unique(dat))!=2) {
    cat("Warning: ", AIM, " is not bi-allelic. Present alleles are ")
    for(i in unique(dat)) cat(i, " ")
    cat("\n")
#    used.AIMs = used.AIMs[used.AIMs != AIM] # determine which AIMs will be used
  }  
}

used.AIMs = used.AIMs[used.AIMs != "rs1919550"]
used.AIMs = used.AIMs[used.AIMs != "rs2024566"]
#used.AIMs = used.AIMs[used.AIMs != "rs16891982"] # This SNP separates EUR from the rest

# exclude non-bi-allelic SNPs
training.sample = training.sample[,used.AIMs]
test.sample = test.sample[,used.AIMs]

# translate training.sample and test.sample to number of occurrences of reference allele
training.sample.numeric = matrix(NA, ncol= ncol(training.sample), nrow = nrow(training.sample))
rownames(training.sample.numeric) = rownames(training.sample)
colnames(training.sample.numeric) = colnames(training.sample)
test.sample.numeric = matrix(NA, ncol= ncol(test.sample), nrow = nrow(test.sample))
rownames(test.sample.numeric) = rownames(test.sample)
colnames(test.sample.numeric) = colnames(test.sample)

for(AIM in colnames(training.sample)) {
  for(i in 1:nrow(training.sample)) {
    if(!is.na(training.sample[i, AIM])) {
      loc = as.character(training.sample[i, AIM])
      loc = unlist(strsplit(loc, ","))
      loc = unlist(strsplit(loc, ""))
      training.sample.numeric[i, AIM] = sum(loc == ref.alleles[AIM])
    }
  }
  for(i in 1:nrow(test.sample)) {
    if(!is.na(test.sample[i, AIM])) {
      loc = as.character(test.sample[i, AIM])
      loc = unlist(strsplit(loc, ","))
      loc = unlist(strsplit(loc, ""))
      test.sample.numeric[i, AIM] = sum(loc == ref.alleles[AIM])
    }
  }
}

# Compute predictions as in Snipper, based on superpops
training.freqs = getfreqs(training.sample.numeric, training.samLoc)

training.refcounts = get.refallelecounts(training.sample.numeric, training.samLoc)
training.allcounts = get.allallelecounts(training.sample.numeric, training.samLoc)

predictions = matrix(0, ncol = nrow(training.refcounts), nrow = nrow(test.sample))
colnames(predictions) = rownames(training.refcounts)
rownames(predictions) = rownames(test.sample)

for(i in 1:nrow(test.sample)) {
  predictions[i,] = prediction.naiveBayes(test.sample.numeric[i,], training.refcounts, training.allcounts) 
}

# Compute admixture and recent-admixture, based on superpops
ia = get.admixtureproportions.multi(test.sample.numeric, training.freqs, tol=1e-6, verbose = TRUE) 
pia = get.recentadmixtureproportions.multi(test.sample.numeric, training.freqs, tol=1e-6, verbose = TRUE) 
loglik.ia = loglik.pia = test.sample.numeric
for(i in 1:nrow(test.sample)) {
  loglik.ia[i,] = loglik.admixture(test.sample.numeric[i,], training.freqs, ia[i,], sum = FALSE)
  loglik.pia[i,] = loglik.recentadmixture(test.sample.numeric[i,], training.freqs, pia[i,], sum=FALSE)
}
delta.loglik = loglik.pia - loglik.ia

save(training.freqs, predictions, ia, loglik.ia, pia, loglik.pia, delta.loglik, test.sample, training.samLoc, file = "analysis_illumina")


