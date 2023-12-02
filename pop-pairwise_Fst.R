# R script to estimate population-pairwise Fst values,
# along with their p-values and 95% confidence intervals.
# Input file needs to be in VCF format.

# needed packages:
#install.packages("vcfR")
#install.packages("dartR")
#install.packages("StAMPP")
#install.packages("reshape2")
library(vcfR)
library(dartR)
library(StAMPP)
library(reshape2)


# set parameters

# set input file name:
input_file_name = "sicula585x39k_m05R7r6h6DP420_1p0p.vcf"

# set output file prefix
output_file_name = "sicula585x39k_m05R7r6h6DP420_1p0p"

# set how many characters of the sample name belong to the population code
poplength = 3

# set number of iterations to run to calculate the p-value for the HWe test
bootrep = 9999

# OPTIONAL: set the order of populations as you want them to appear in the output tables
# if you want to keep the defualt (alphabetical) order set to NULL
poporder <- NULL


#######################################################################


# import VCF file
myvcf <- read.vcfR(input_file_name)
myvcf

# extract sample names
mysamples <- (dimnames(myvcf@gt)[2][[1]])[-1]

# extract population IDs
mypops <- substr(mysamples, 1, poplength)

# extract a total number of populations
maxpop <- length(levels(as.factor(mypops)))

# transform VCF to GENIND format
mygenind <- vcfR2genind(myvcf)

# add population IDs to GENIND file
mygenind@pop <- as.factor(mypops)

# transform GENIND format to GENLIGHT format
mygenlight <- gi2gl(mygenind)

# run pairwise Fst calculation
pairFst <- stamppFst(mygenlight, nboots=bootrep, percent=95, nclusters=1)

# extract Fst and p-value matrix 
fsts <- pairFst$Fsts
pvals <- pairFst$Pvalues

# extract confidence intervals
cis <- pairFst[["Bootstraps"]][,c("Population1", "Population2", "Lower bound CI limit", "Upper bound CI limit")]
cis_lower <- as.data.frame(reshape2::dcast(cis, Population1~Population2 , value.var='Lower bound CI limit', fill=NA))
cis_upper <- as.data.frame(reshape2::dcast(cis, Population1~Population2 , value.var='Upper bound CI limit', fill=NA))

# edit the confidence intervals tables
colnames(cis_lower)[1] <- cis_lower[1,1]
colnames(cis_upper)[1] <- cis_upper[1,1]
row.names(cis_lower) <- cis_lower[,1]
row.names(cis_upper) <- cis_upper[,1]
cis_lower[,1] <- NA
cis_upper[,1] <- NA
cis_lower[maxpop,] <- NA
cis_upper[maxpop,] <- NA
rownames(cis_lower)[maxpop] <- colnames(cis_lower)[maxpop]
rownames(cis_upper)[maxpop] <- colnames(cis_upper)[maxpop]
cis_lower <- as.matrix(cis_lower)
cis_upper <- as.matrix(cis_upper)

# make the matrices symmetrical
fsts[upper.tri(fsts)]  <- t(fsts)[upper.tri(fsts)]
pvals[upper.tri(pvals)]  <- t(pvals)[upper.tri(pvals)]
cis_lower[lower.tri(cis_lower)]  <- t(cis_lower)[lower.tri(cis_lower)]
cis_upper[lower.tri(cis_upper)]  <- t(cis_upper)[lower.tri(cis_upper)]

# reorder the matrices if so specified
if (is.null(poporder) == FALSE) {
  
  fsts <- fsts[poporder,poporder]
  pvals <- pvals[poporder,poporder]
  cis_lower <- cis_lower[poporder,poporder]
  cis_upper <- cis_upper[poporder,poporder]
}

# create an output file
fstfile = paste("pairFst_", output_file_name, ".txt", sep="")

# export results
cat("Fsts:", sep="\n", file=fstfile, append=FALSE)
write.table (format (fsts, digits=5), file=fstfile, sep="\t", append=TRUE, quote=FALSE)
cat("\n", "p-values: ", "\n", sep="",file=fstfile, append=TRUE)
write.table (format (pvals, digits=5), file=fstfile, sep="\t", append=TRUE, quote=FALSE)
cat("\n", "lower 95% CI: ", "\n", sep="",file=fstfile, append=TRUE)
write.table (format (cis_lower, digits=5), file=fstfile, sep="\t", append=TRUE, quote=FALSE)
cat("\n", "upper 95% CI: ", "\n", sep="",file=fstfile, append=TRUE)
write.table (format (cis_upper, digits=5), file=fstfile, sep="\t", append=TRUE, quote=FALSE)

