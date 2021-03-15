# recent-admixture
This repository is accompanying the paper
Inference of recent admixture using genotype data by Peter Pfaffelhuber, Elisabeth Huss, Franz Baumdicker, Jana Naue, Sabine Lutz-Bonengel, Fabian Staubach.

The main functions are found in code/tools.r. The variables are as follows:
ninds: number of individuals
nss: number of segregating sites
npops: number of superpopulations
filenameDataVCFGZ is a file of format vcf.gz
inds: ninds x 3 matrix, where the first column are the column names from the vcf.gz file which we want to read; the second column is the population, the third the superpopulation
samLoc is a vector of length ninds, which stores superpopulations; names are the individuals
sample has the genotypes (0 and 2 for homozygotes, 1 for heterozygote, NA for not available); the SNPs are in the columns, individuals in the rows.
refallelecounts.training: npops x nss matrix, where entry [i,j] counts the occurrences of the derived allele in population i for SNP j.
allallelecounts.training: npops x nss matrix, where entry [i,j] counts the occurrences of available alleles in population i for SNP j.
freqs: npops x nss matrix, where entry [i,j] is the allele frequency of the derived allele in population i for SNP j.
admixture.proportions: vector of length npops, where entry i is the proportion of the genome which comes from population i. 
recentadmixture.proportions: vector of length 2*npops, where entries i and npops+i are the proportions of the genomes of the parents which come from population i. 


getData(filenameDataVCFGZ, inds):
The result is a list with variables sample (the data) and samLoc (the locations)

prediction.naiveBayes(sample.test, refallelecounts.training, allallelecounts.training) 
Classification through a naive Bayes classifier. Result is a matrix of size ninds x npops, where entry [i,j] gives the class probability of individual i for population j. This function needs
log.prediction.naiveBayes(sample.test, refallelecounts.training, allallelecounts.training),
get.refallelecounts = function(sample, samLoc)
get.allallelecounts = function(sample, samLoc).

getfreqs = function(sample, samLoc)
Computes allele frequencies in all populations. Result is a matrix of size npops x nss

get.admixtureproportions.multi(sample, freqs, tol=1e-6, verbose = FALSE) 
For the admixture model, this function computes individual admixture (IA) based on the allelic frequencies in all populations. Computation is an iteration. It stops when the iteration does not produce a difference of more than tol. verbose determines the output. Here, sample is a matrix. Result is a nrow(sample) x npops matrix, where entry [i,j] is the proportion of the genome of individual i which comes from population j.

get.admixtureproportions(sample, freqs, tol=1e-6, verbose = FALSE)
Same as ...multi, but sample is a vector, i.e. only a single individual is studied. Result is a vector of length npops.

get.recentadmixtureproportions.multi(sample, freqs, tol=1e-6, verbose = FALSE)
Same as ...admixtureproportions..., but for the recent-admixture model (which assumes different individual admixture for both parents, i.e. parental individual admixutre, PIA). Result is a nrow(sample) x (2*npops) matrix, where entries [i,j] and [i, npops+j] are the proportions of the genome of the parents of individual i which come from population j.

get.recentadmixtureproportions(sample, freqs, tol=1e-6, verbose = FALSE)
Same as ...multi, but sample is a vector, i.e. only a single individual is studied. Result is a vector of length 2*npops.

fun3(ia, freqs, loc.sample)
This function is called iteratively from get.recentadmixtureproportions.

fun2(ia, freqs, loc.sample) 
This function is called iteratively from get.admixtureproportions.

get.loglik.admixture(loc.sample, freqs, admixture.proportions, sum=TRUE)
Log-likelihood function for the admixture model. If sum=TRUE, a single value is reported. If sum=FALSE, the contributions of the single SNPs to the log-likelihood are given. 

get.loglik.recentadmixture(loc.sample, freqs, recentadmixture.proportions, sum = TRUE) 
Same as loglik.admixture, but for the recent-admixture model.

rdirichlet(n, classes)
Simulates a Dirichlet(1,...,1)-distributed random variable, where classes is the number of entries and n is the number of simulated vectors. This is used as starting point to get.admixtureproportions and  get.recentadmixtureproportions.

simulate.p.value(loc.sample, freqs, nsam=100, tol=1e-06, verbose = TRUE) 
For simulation of p-values for the trace loc.sample and the allele frequencies of the reference database freqs, this function simulates nsam individuals with the same admixtureproportions as loc.sample, and reports back the frequency of delta-loglikelihoods above the trace.


In data/1000G, all relevant data and analysis scripts from the 1000 genomes dataset are stored.
1000G_AIMsetEUROFORGEN.vcf.gz: Data for 121 SNPs from the EUROFORGEN AIMset in the 1000 genomes dataset. 
1000G_AIMsetKidd.vcf.gz: Same for the Kidd AIMset
1000G_SampleListWithLocations.txt: File containing sampling locations in the 1000 genomes dataset
analysis_EUROFORGENE.r: Script for the analysis of the 1000 genomes data when using the EUROFORGEN AIMset.
analysis_Kidd.r: Same for the Kidd AIMset.

In data/self, the scripts for analysing the self-collected individuals are stored, as well as the reference datasets.
Note that the required file SNPstwoVPs.csv which contains genotype information is not contained in the repository. Please ask the authors when you are interested in this data.

If you want to analyse your own data, you can use the files in data/self as well. Then, you have to provide a csv-file with the same structure as SNPstwoVPs.csv (ideally with tabs as delimiter) and change line 23 in analysis_iullumina.r in order to read this file.

The R-scripts do not produce figures. Results are variables like (for the self-collected individuals):
training.freqs: Allele frequencies in the training set
predictions: For the naive Bayes classifier, prediction probabilities are given.
ia: Individual Admixture for all samples
loglik.ia: Log-likelihood of the admixture modek for all samples, broken down for all SNPs
pia: Parental Individual Admixture for all samples (each sample is one row, the two parents are reported consecutively)
loglik.pia: Log-likelihood of the recent-admixture modek for all samples, broken down for all SNPs
delta.loglik: Difference of rowSums(loglik.pia) and rowSums(loglik.pia)
test.sample: The data which is analysed
training.samLoc: Sampling locations for individuals in the training set

