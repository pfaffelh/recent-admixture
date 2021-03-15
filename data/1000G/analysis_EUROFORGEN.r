require("vcfR")
source("../../code/tools.r")
# read reference data from Snipper input file
nsam = 2000 # samples for simulating p-values
require("parallel")
options(mc.cores = 3)

##################

a = read.table("1000G_SampleListWithLocations.txt", header = FALSE)
inds = a[,1:3]
dat = getData("1000G_AIMsetEUROFORGEN.vcf.gz", inds)
test.sample.numeric = dat$sample
test.samLoc = dat$samLoc
test.samLocFine = dat$samLocFine
names(test.samLoc) = rownames(test.sample.numeric)
names(test.samLocFine) = rownames(test.sample.numeric)

training.indsUsed = !is.element(test.samLocFine, c("ASW", "ACB", "MXL", "PUR", "CLM"))
training.sample.numeric = test.sample.numeric[training.indsUsed,]
training.samLoc = test.samLoc[training.indsUsed]
training.samLocFine = test.samLocFine[training.indsUsed]
names(training.samLoc) = rownames(training.sample.numeric)
names(training.samLocFine) = rownames(training.sample.numeric)

# select the three individuals treated in the paper
test.sample.numeric = test.sample.numeric[c("NA19719", "NA19720", "NA20278"),]
test.samLoc = test.samLoc[c("NA19719", "NA19720", "NA20278")]
test.samLocFine = test.samLocFine[c("NA19719", "NA19720", "NA20278")]

# Compute predictions as in Snipper, based on superpops
training.freqs = getfreqs(training.sample.numeric, training.samLoc)

training.refcounts = get.refallelecounts(training.sample.numeric, training.samLoc)
training.allcounts = get.allallelecounts(training.sample.numeric, training.samLoc)

predictions = matrix(0, ncol = nrow(training.refcounts), nrow = nrow(test.sample.numeric))
colnames(predictions) = rownames(training.refcounts)
rownames(predictions) = rownames(test.sample.numeric)

for(i in 1:nrow(test.sample.numeric)) {
  predictions[i,] = prediction.naiveBayes(test.sample.numeric[i,], training.refcounts, training.allcounts) 
}

# Compute admixture, recent-admixture and p-value, based on superpops

ia = get.admixtureproportions.multi(test.sample.numeric, training.freqs, tol=1e-6, verbose = TRUE) 
pia = get.recentadmixtureproportions.multi(test.sample.numeric, training.freqs, tol=1e-6, verbose = TRUE) 
loglik.ia = loglik.pia = test.sample.numeric
p = test.sample.numeric[,1]

res = mclapply(1:nrow(test.sample.numeric), function(i){
    cat("\n\nIndividual ", i, "; Calculate p-value for ", rownames(test.sample.numeric)[i], " from ", as.vector(test.samLoc[i]), ".\n", sep="")
    res = simulate.p.value(test.sample.numeric[i,], training.freqs, nsam, tol=1e-06, verbose = TRUE)
    cat("p-value = ", res$p, ".\n", sep="")
    res
  })


for(i in 1:nrow(test.sample.numeric)) {
  loglik.ia[i,] = get.loglik.admixture(test.sample.numeric[i,], training.freqs, ia[i,], sum = FALSE)
  loglik.pia[i,] = get.loglik.recentadmixture(test.sample.numeric[i,], training.freqs, pia[i,], sum=FALSE)
  p[i] = res[[i]]$p
}
delta.loglik = loglik.pia - loglik.ia

save(training.freqs, predictions, ia, loglik.ia, pia, loglik.pia, delta.loglik, test.sample.numeric, training.samLoc, file = "analysis_EUROFORGEN")

for(i in 1:nrow(test.sample.numeric)) {
  cat("\n\nAnalysis for individual ", names(p)[i], ".\n", sep="")
  cat("ia ", round(as.matrix(t(ia[i,])),3), "\n")
  cat("pia ", round(as.matrix(t(pia[i,])),3), "\n")
  cat("Heterozygous sites ", sum(test.sample.numeric[i,]==1), ".\n", sep="")
  cat("Contribution to delta ell ", round(sum(delta.loglik[i, test.sample.numeric[i,]==1]),3), ".\n\n", sep="")
  cat("Homozygous sites ", sum(test.sample.numeric[i,]!=1), ".\n", sep="")
  cat("Contribution to delta ell ", round(sum(delta.loglik[i, test.sample.numeric[i,]!=1]),3), ".\n\n", sep="")
  cat("Sum of sites ", sum(test.sample.numeric[i,]<10), ".\n", sep="")
  cat("Contribution to delta ell ", round(sum(delta.loglik[i,]),3), ".\n\n", sep="")
  cat("p-value ", 2*p[i], ".\n\n", sep="")
}

