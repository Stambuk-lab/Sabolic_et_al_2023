# R script to estimate genomic diversity indices,
# along with their p-values and 95% confidence intervals.
# Input file needs to be in GENPOP format.

# needed packages:
#install.packages("diveRsity")
#devtools::install_github('wrengels/HWxtest', subdir='pkg')
library(diveRsity)


# set parameters

# set input file name:
genepopfile = "sicula585x39k_m05R7r6h6DP420_1p0p.vcf"

# set output file prefix
output_file_name = "sicula585x39k_m05R7r6h6DP420_1p0p"

# set the number of bootstrap replicates to generate confidence intervals for Fis and Ar
bootci = 999

# set the number of bootstrap iterations to generate MCMC for HWe test
bootrep = 9999



#######################################################################


# get genomic diversity indices
gen_div <- basicStats(infile = genepopfile, outfile = NULL, 
                      fis_ci=TRUE, ar_ci=TRUE, fis_boots = bootci, ar_boots = bootci, 
                      mc_reps = bootrep, rarefaction=FALSE, ar_alpha=0.05, fis_alpha=0.05)

# extract main results 
totals <- lapply(gen_div$main_tab, function(x) (x)$overall)

# edit population IDs
popnames <- gsub('[[:digit:]]+', '', names(totals))
names(totals) <- popnames

# transform the output to a dataframe
totals.df <- as.data.frame(totals)

# add indies IDs
indexnames <- rownames(gen_div$main_tab[[1]])
rownames(totals.df) <- indexnames

# transpose the dataframe
indexes <- as.data.frame(t(totals.df))

# extract result table
div_stats = paste("genomic_diversity_stats_", output_file_name, ".txt", sep="")
cat("Diversity indexes (diveRsity): ", "", sep="\n", file=div_stats, append=FALSE)
write.table (format(indexes, digits=4), file=div_stats, sep="\t", append=TRUE, quote=FALSE)

