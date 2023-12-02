# R script to find and output a list of loci not in HWe in a given dataset,
# both globally (across the entire dataset)
# and across specified n% of populations in the dataset.
# Input file needs to be in VCF format.

# needed packages:
#install.packages("vcfR")
#install.packages("adegenet")
#install.packages("pegas")
library(vcfR)
library (adegenet)
library(pegas)


# set parameters

# set input file name:
input_file_name = "sicula585x39k_m05R7r6h6DP420_1p0p.vcf"

# set output file prefix
output_file_name = "sicula585x39k_m05R7r6h6DP420_1p0p"

# set how many characters of the sample name belong to the population code
poplength = 3

# set number of iterations to run to calculate the p-value for the HWe test
iterations = 999

# set threshold (critical) value below which a pvalue is considered significant
refval = 0.05

# set proportion of populations in which a locus must be out of HWE in order for it to be deleted
ref_outhw = 0.6


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

# extract a list of loci/SNP names
locilist <- myvcf@fix[,"ID"]

# extract total number of loci/SNPs in the input file
locinum <- length(locilist)

# transform VCF to GENIND format
mygenind <- vcfR2genind(myvcf)

# add popnames to GENIND file
mygenind@pop <- as.factor(mypops)

# split GENIND file by population
mygenind_popdata <- seppop(mygenind)

# build a output dataframe with locinames on the rows
basetable <- as.data.frame(locilist)

# create empty output lists
pophw <- list()
signifpval <- list()
logfile = paste0(output_file_name, "_n", ref_outhw, "_popHWe_log.txt")

# run HWE test
for (population in mygenind_popdata) { 
  
  # print which population the analysis is run on
  popname <- levels(population@pop)
  cat("Running HWe test for:", popname, "\n")

  # run HWe test
  pophw[[popname]] <- as.data.frame(hw.test(population, B=iterations))
  
  # extract SNPs out of HWe
  signifpval[[popname]] <- rownames(pophw[[popname]])[which(pophw[[popname]][["Pr.exact"]] < refval)]

  # print statistics for each population
  signum <- length(signifpval[[popname]])
  proportionsig = signum / locinum
  percentsig = round(proportionsig*100, digits=2)
  cat(popname, ": ", signum, " out of ", locinum, " loci were significantly out of HWe (", percentsig, "% p-values were below ", refval, ")\n\n", sep="")
  
  # extract statistics to log file
  cat(popname, ": ", signum, " out of ", locinum, " loci were significantly out of HWe (", percentsig, "% p-values were below ", refval, ")\n\n", sep="", file=logfile, append=TRUE)
  
  # subset and extract p-values
  infoHWtable <- data.frame(pophw[[popname]][["Pr.exact"]], rownames(pophw[[popname]]))
  names(infoHWtable )<- c(popname, "locilist")

  # merge to output dataframe
  basetable <- merge(basetable, infoHWtable, by="locilist")
}

# extract only p-value columns from the basetable dataframe
lociinfo <- basetable[,c(2:(maxpop+1))]

# add columns to basetable to hold the number of pops in which the loci is out of hwe
basetable$pophw <- NA
basetable$outhw <- NA

# create an empty output list that will contain SNP ids of loci that are out of HWe in n% of populations
filterhw_loci_list = c()

# loop and add the number of populations in which each loci is out of hwe
for (loci in 1:length(locilist)) {
  
  # extract a number of pops in whic locus is out of HWe
  outnum <- sum(lociinfo[loci,] < refval)
 
  # extract a proportion of pops in whic locus is out of HWe in regards to the total number of pops
  hwi_rate <- (outnum/maxpop)
  
  # add to the table
  basetable[loci, "pophw"] <- outnum
  basetable[loci, "outhw"] <- hwi_rate
  
  #if proportion is too high, add locus name to delete list
  locusname <- as.character(basetable[loci,1])
  if (hwi_rate >= ref_outhw) { 
    filterhw_loci_list <- c(filterhw_loci_list, locusname)
  }
}

# extract statistics
filterlength <- length(filterhw_loci_list)
finalnum = locinum - filterlength

cat("The \"worst\" locus is out of HWe in", max(basetable$outhw)*maxpop, "out of", maxpop, "populations.", "\n")
cat("That is a frequency of", max(basetable$outhw), "across all populations; only loci above", ref_outhw, "were considered out of HWe.", "\n")
cat("In total", filterlength, "out of", locinum, "loci are out of HWe.", "\n")
cat("If you filter them out, there will be", finalnum, "loci left.")

# export the basetable dataframe
write.table(basetable, file=paste0(output_file_name, "_n", ref_outhw, "_popHWe_stats.txt"), 
            sep="\t", col.names=T, row.names=F, quote=F)

# export a list of loci out of HWe to be filtered out
write.table(filterhw_loci_list, file=paste0(output_file_name, "_n", ref_outhw, "_popHWe.txt"), 
            sep="\t", col.names=F, row.names=F, quote=F)

# run global HWe test
cat("Running global HWe test of the entire dataset...")
globalhw <- as.data.frame(hw.test(mygenind, B=iterations))

# extract SNPs out of HWe globally
signifpval_global <- rownames(globalhw[which(globalhw$Pr.exact <refval),])

# extract statistics
filterlength <- length(signifpval_global)
finalnum = locinum - filterlength
cat("In total", filterlength, "out of", locinum, "loci are out of HWe globally.", "\n")
cat("If you filter them out, there will be", finalnum, "loci left.")

# export the the main result dataframe
write.table(globalhw, file=paste0(output_file_name, "_globalHWe_stats.txt"), 
            sep="\t", col.names=T, row.names=T, quote=F)

# export a list of loci out of HWe to be filtered out
write.table(signifpval_global, file=paste0(output_file_name, "_globalHWe.txt"), 
            sep="\t", col.names=F, row.names=F, quote=F)

