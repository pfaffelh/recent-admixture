require("vcfR")
require("parallel")

# inds is nsam x 3 matrix, where the first column are the column names from the vcf.gz file which we want to read;
# the second column is the population, the third the superpopulation

getData<-function(filenameDataVCFGZ, inds) {
  trans = c(0,1,1,2, NA)
  names(trans)=c("0|0", "0|1", "1|0", "1|1", "NA")
  samLocFine = inds[,2]
  samLoc = inds[,3]
  
  a = read.vcfR(filenameDataVCFGZ, skip = 0, nrows = -1)
  genotypes = a@gt[,-1]
  row.names(genotypes) = a@fix[,3]
  sample = matrix(0, ncol = nrow(genotypes), nrow = nrow(inds))
  rownames(sample) = as.character(inds[,1])
  colnames(sample) = as.character(a@fix[,3])
  for(k in 1:nrow(inds)) {
    sample[k,] = as.numeric(trans[genotypes[,as.character(inds[k,1])]])
  }
  res = list(sample = sample, samLoc = samLoc, samLocFine = samLocFine)
  res
}

# Here is a naive Bayesian classifier, similar to Snipper
# sample.test is a matrix
# We use prediction.naiveBayes(...) = exp(log.prediction.naiveBayes(...)) for numerical stability
# sample.test is a vector

log.prediction.naiveBayes<-function(sample.test, refallelecounts.training, allallelecounts.training) {
  refallelecounts.training = refallelecounts.training[,!is.na(sample.test)]
  allallelecounts.training = allallelecounts.training[,!is.na(sample.test)]
  sample.test = sample.test[!is.na(sample.test)]

  demes = nrow(allallelecounts.training)
  res = 1:demes
  names(res) = rownames(allallelecounts.training)

  for(k in 1:demes) {
    res[k] = sum(log((1 + refallelecounts.training[k,])) * sample.test + log((1 + allallelecounts.training[k,] - refallelecounts.training[k,])) * (2 -t(sample.test)) - log(1 + allallelecounts.training[k,]) * 2)
  } 
  res 
}

# See log.prediction.naiveBayes;

prediction.naiveBayes<-function(sample.test, refallelecounts.training, allallelecounts.training) {
  res = log.prediction.naiveBayes(sample.test, refallelecounts.training, allallelecounts.training) 
  res = exp(res - max(res))
  res / sum(res)
}

# functions to get the frequencies of SNPs per deme

get.refallelecounts = function(sample, samLoc){
  t(simplify2array(by(sample,samLoc,colSums,na.rm=TRUE)))
}

get.allallelecounts = function(sample, samLoc){
  sample[!is.na(sample)] = 2
  t(simplify2array(by(sample,samLoc,colSums,na.rm=TRUE))) 
}

# compute frequencies from sample along samLoc
getfreqs = function(sample, samLoc){
  t(simplify2array(by(sample,samLoc,colMeans,na.rm=TRUE))) / 2
}

# sample has nss columns and ninds rows; entries are 0, 1, 2
# freqs has npops rows and nss columns
# countsAncestral is a vector of length npops
# value is a matrix admixture.proportions with nsam.mixed rows and npops columns;
# admixture.proportion[i,k] gives the proportion of the genome of individual i coming from population k

get.admixtureproportions.multi<-function(sample, freqs, tol=1e-6, verbose = FALSE) {
  res = NULL
  for(i in 1:nrow(sample)) {
    if(verbose) cat("Method: admixture \t Name: ", rownames(sample)[i])
    res = rbind(res, get.admixtureproportions(sample[i,], freqs, tol=tol, verbose=verbose))
  }
  colnames(res) = rownames(freqs)
  rownames(res) = rownames(sample)
  res
}

# Here, sample is a vector

get.admixtureproportions<-function(sample, freqs, tol=1e-6, verbose = FALSE) {
  freqs = freqs[,!is.na(sample)]
  sample = sample[!is.na(sample)]
  
  npops = nrow(freqs)
  nss = ncol(freqs)
  
  # initialize admixture proportions
  res = as.vector(rdirichlet(1,npops))
  #  res = rep(1/npops, npops)
  names(res) = rownames(freqs)

  err = 1
  j=1
  while(err>tol) {
    j = j+1
    loc = fun2(res, freqs, sample)
    err = sum(abs(res - loc))
    res = loc
  }
  if(verbose) cat("\t Iterations: ", j, "\n")
  res
}

get.recentadmixtureproportions.multi<-function(sample, freqs, tol=1e-6, verbose = FALSE) {
  res = NULL
  for(i in 1:nrow(sample)) {
    if(verbose) cat("Method: recent-admixture \t Name: ", rownames(sample)[i])
    res = rbind(res, get.recentadmixtureproportions(sample[i,], freqs, tol=tol, verbose=verbose))
  }
  colnames(res) = c(rownames(freqs), rownames(freqs))
  rownames(res) = rownames(sample)
  res
}

get.recentadmixtureproportions<-function(sample, freqs, tol=1e-6, verbose = FALSE) {
  freqs = freqs[,!is.na(sample)]
  sample = sample[!is.na(sample)]
  npops = nrow(freqs)
  nss = ncol(freqs)
  
  # initialize admixture proportions, res has two rows for the two parents
#  a = get.admixtureproportions(sample, freqs, tol, verbose = FALSE)
#  a = rep(1/npops, npops)
#  a = rep(1/npops, npops)
#  res = rbind(a,a) # matrix(1/npops, nrow=2, ncol=npops)
  res = rdirichlet(2,npops)
  colnames(res) = rownames(freqs)

  err = 1
  j=1
  while(err>tol) {
    j = j+1
    loc = fun3(res, freqs, sample)
    err = sum(abs(res - loc))
    res = loc
    }
  if(verbose) cat("\t Iterations: ", j, "\n")
  c(res[1,], res[2,])
}

fun3<-function(ia, freqs, loc.sample) {
  npops = ncol(ia)
  nss = length(loc.sample)  
  E = matrix(0, nrow=npops, ncol = nss)

  loc = ia %*% freqs # has nss cols and two rows
  loc[loc==0] = 1e-16
  loc[loc==1] = 1-1e-16
  for(k in 1:npops) {
    E[k,] = (loc.sample==2)*freqs[k,]/loc[2,] + (loc.sample==0)*(1-freqs[k,])/(1-loc[2,]) + (loc.sample==1)*(freqs[k,]*(1-loc[1,]) + (1-freqs[k,])*loc[1,])/(loc[1,]*(1-loc[2,]) + loc[2,]*(1-loc[1,]))
  }
  ia[2,] = rowSums(E)/(nss) * ia[2,]
  ia[2,] = ia[2,]/(sum(ia[2,]))

  for(k in 1:npops) {
    E[k,] = (loc.sample==2)*freqs[k,]/loc[1,] + (loc.sample==0)*(1-freqs[k,])/(1-loc[1,]) + (loc.sample==1)*(freqs[k,]*(1-loc[2,]) + (1-freqs[k,])*loc[2,])/(loc[2,]*(1-loc[1,]) + loc[1,]*(1-loc[2,]))
  }
  ia[1,] = rowSums(E)/(nss) * ia[1,]
  ia[1,] = ia[1,]/(sum(ia[1,]))

  ia
}

fun2<-function(ia, freqs, loc.sample) {
  npops = length(ia)
  nss = length(loc.sample)  
  E = matrix(0, nrow=npops, ncol = nss)

  # admixture.proportions = admixture.proportions / sum(admixture.proportions)

  loc = ia %*% freqs # has nss cols
  loc[loc==0] = 1e-16
  loc[loc==1] = 1-1e-16
  for(k in 1:npops) {
    E[k,] = (loc.sample * freqs[k,] / loc + (2 - loc.sample) * (1-freqs[k,])/(1 - loc))
  }
  res = rowSums(E)/(2*nss) * ia
  res/(sum(res))
}

get.loglik.admixture<-function(loc.sample, freqs, admixture.proportions, sum=TRUE) {
  loc = t(freqs) %*% admixture.proportions
  if(sum) res = sum(log(choose(2, loc.sample) * loc^loc.sample * (1-loc)^(2-loc.sample)))
  else res = log(choose(2, loc.sample) * loc^loc.sample * (1-loc)^(2-loc.sample))
  res
}

get.loglik.recentadmixture<-function(loc.sample, freqs, recentadmixture.proportions, sum = TRUE) {
  if(is.vector(recentadmixture.proportions)) {
    n = length(recentadmixture.proportions)
    recentadmixture.proportions = rbind(recentadmixture.proportions[1:(n/2)], recentadmixture.proportions[(n/2+1):n])
  }  
  loc1 = t(freqs) %*% recentadmixture.proportions[1,]
  loc2 = t(freqs) %*% recentadmixture.proportions[2,]
  if(sum) res = sum(log((loc.sample==2) * loc1 * loc2 + (loc.sample==1) * (loc1 * (1-loc2) + (1-loc1)*loc2) + (loc.sample==0) * (1-loc1) * (1-loc2)))
  else res = log((loc.sample==2) * loc1 * loc2 + (loc.sample==1) * (loc1 * (1-loc2) + (1-loc1)*loc2) + (loc.sample==0) * (1-loc1) * (1-loc2))
  res
}

rdirichlet<-function(n, classes){
  res = matrix(rexp(classes*n), ncol=classes, nrow=n)
  res / rowSums(res)
}

simulate.p.value<-function(loc.sample, freqs, nsam=100, tol=1e-06, verbose = TRUE) {
  loc.ia<-get.admixtureproportions(loc.sample, freqs, tol, verbose = FALSE)
  loc.pia<-get.recentadmixtureproportions(loc.sample, freqs, tol, verbose = FALSE)
  loc.loglik.recentadmixture = get.loglik.recentadmixture(loc.sample, freqs, loc.pia)
  loc.loglik.admixture = get.loglik.admixture(loc.sample, freqs, loc.ia)
  loc.delta.ell = loc.loglik.recentadmixture - loc.loglik.admixture

  demes = nrow(freqs)
  if(verbose) {
    cat("Simulating p-value for delta.ell =", loc.delta.ell, "and pia =", round(loc.pia,3), ".\n", sep=" ")
  } 
  # simulate nsam individuals with the same ia; store result in sample
  sample = matrix(0, ncol=ncol(freqs), nrow = nsam)
  rownames(sample) = 1:nsam
  colnames(sample) = colnames(freqs)
  for(j in 1:ncol(freqs)) {
    p = loc.ia %*% freqs[,j]
    sample[,j] = rbinom(nsam, 2, p)
  }
  ia =  get.admixtureproportions.multi(sample, freqs, tol=tol, verbose=FALSE)
  pia = get.recentadmixtureproportions.multi(sample, freqs, tol=tol, verbose=FALSE)
  delta.ell = sapply(1:nrow(sample), function(x) get.loglik.recentadmixture(sample[x,], freqs, loc.pia) - get.loglik.admixture(sample[x,],
     freqs, loc.ia))
  res = list(ia = loc.ia, pia = loc.pia, loglik.recentadmixture = loc.loglik.recentadmixture, loglik.admixture = loc.loglik.admixture,
    delta.ell = loc.delta.ell, p = mean(delta.ell>loc.delta.ell))
  res
}

