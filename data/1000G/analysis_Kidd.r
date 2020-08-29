require("vcfR")
source("../../code/tools.r")
# read reference data from Snipper input file

##################

a = read.table("1000G_SampleListWithLocations.txt", header = FALSE)
inds = a[,1:3]
dat = getData("1000G_AIMsetKidd.vcf.gz", inds)
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

# Compute admixture and recent-admixture, based on superpops
ia = get.admixtureproportions.multi(test.sample.numeric, training.freqs, tol=1e-6, verbose = TRUE) 
pia = get.recentadmixtureproportions.multi(test.sample.numeric, training.freqs, tol=1e-6, verbose = TRUE) 
loglik.ia = loglik.pia = test.sample.numeric
for(i in 1:nrow(test.sample.numeric)) {
  loglik.ia[i,] = loglik.admixture(test.sample.numeric[i,], training.freqs, ia[i,], sum = FALSE)
  loglik.pia[i,] = loglik.recentadmixture(test.sample.numeric[i,], training.freqs, pia[i,], sum=FALSE)
}
delta.loglik = loglik.pia - loglik.ia

save(training.freqs, predictions, ia, loglik.ia, pia, loglik.pia, delta.loglik, test.sample.numeric, training.samLoc, file = "analysis_Kidd")


