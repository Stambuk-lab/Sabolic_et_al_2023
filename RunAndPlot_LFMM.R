# R script for for running LFMM analysis using lfmm() package,
# and extracting and plotting the results.
# It will:
# run lfmm analysis using ridge penalty method
# and a range of modified GIF values,
# extract significant SNPs after Benjamini-Hochberg correction of p-values using specified FDR threshold,
# and plot main results for each variable.

# Genotype file should be in VCF format.

# The matrix with additional data to give to model (i.e. phenotype or environmental data)
# needs to look like this (no missing values allowed!):
# individual data -> grouping variables -> additional variables, i.e.:
# INDIVIDUALS  non-inf1  non-inf2  ...  non-infn   add1  add2  ...  addn   
#  1	  	     04        M         ...  1          3.2   2.7   ...  4.2  
#  2           07        P         ...  2          5.7   5.8   ...  3.5  
#  3		       02	       E         ...  2          0.8   1.3   ...  8.4 
# ...etc

# needed packages:
#Rtools may be needed
#install.packages("devtools")
#devtools::install_github("bcm-uga/lfmm")
#install.packages("adegenet")
#BiocManager::install("qvalue")
#install.packages("dartR")
#install.packages("venn")
#install.packages("vcfR")
library(devtools)
library(lfmm)
library(adegenet)
library(qvalue)
library(dartR)
library(venn)
library(vcfR)


# set parameters:
outputname = "PKPM_males_K2"                              # set prefix for all output file names (no extensions)

genotype_file = "NPK-NPM82x10k_m05R7r6h6DP420_1p0p.vcf"      # name of the structure file containing genotype data (with .str extension)
subset_ind = TRUE                                            # do you want to subset the individuals from the genotype file
keep_ind = "popmap_male53x10k_NPK-NPM82x10k.txt"             # list of individuals to keep in the genotype file

variable_file = "malePKPM2016_residuals.txt"                       # filename of the matrix with additional data to give to model (should be tab separated .txt)
groupingv = 3                   # set the number non-informative variables
varableno = 16                  # set the number of additional measured variables (phenotypic) in variable_file

k_value = 2                 # K value, i.e. latent factor value to be used in the analysis
fdr = 0.05                  # FDR (q-value) cutoff for significant loci

mod.gif = c(0.5,1.2)            # set the first and last number of a range of GIF values you want to test
by.gif = 0.01                   # set the interval for a range of GIF values you want to test


############## NO NEED TO SET ANYTHING UNDER THIS LINE ###############

###### import input files ######

# import file with additional data to give to model
mydata <- read.table(variable_file, header = T, sep="\t") 
str(mydata)
colnames(mydata)[1] <- "INDIVIDUALS"

# import VCF file
myvcf <- read.vcfR(genotype_file)
myvcf

# transform VCF to GENLIGHT format
mygenlight <- vcfR2genlight(myvcf, n.cores = 1)
mygenlight

# subset genotype file if so specified
if (subset_ind == TRUE) {
    
    # import data on individuals to subset
    my_keep_ind <- read.table(keep_ind)[,1]
    
    # find individuals to exclude
    my_exclude_ind <-  mygenlight@ind.names[!mygenlight@ind.names %in% my_keep_ind]

    # remove female or male individuals
    mygenlight <- gl.drop.ind(mygenlight, my_exclude_ind, mono.rm=TRUE)
    mygenlight
}

coded_genotypes <- as.data.frame(mygenlight)


###### merge data ######

# edit genotype file so it can be merged
names <- rownames(coded_genotypes)
genonames <- colnames(coded_genotypes)
rownames(coded_genotypes) <- NULL
coded_genotypes <- cbind(names,coded_genotypes)
colnames(coded_genotypes)[1] <- "INDIVIDUALS"

# combine datasets
combdata <- merge(mydata, coded_genotypes, by="INDIVIDUALS")


###### extract info about variables #######
id = 1

firstgroup = id+1
lastgroup = id+groupingv
colgroup = c(firstgroup:lastgroup)
groupnames = colnames(combdata[colgroup])

firstpredict = id+groupingv+1
lastpredict = id+groupingv+varableno
predictors = combdata[,firstpredict:lastpredict]
colpred = c(firstpredict:lastpredict)
prednames = colnames(combdata[colpred])

responses <- combdata[,genonames]

# scale predictors
predictors <- as.data.frame(scale(predictors))
colnames(predictors) = prednames


###### run the analysis ######

# fit a LFMM, i.e, compute B, U, V estimates
mod.lfmm <- list()

for (i in prednames) {

    mod.lfmm[[i]] <- lfmm_ridge(Y = responses, X = predictors[i], K = k_value, lambda = 1e-05,
                                algorithm = "analytical", it.max = itno, relative.err.min = relerr)

}

# perform association testing using the fitted model
pv <- list()

for (i in 1:length(mod.lfmm)) {
    
    varname <- names(mod.lfmm[i])
    pv[[varname]] <- lfmm_test(Y = responses, X = predictors[varname], lfmm = mod.lfmm, calibrate = "gif")
}

# extract obtained genomic inflation factors
filename = paste0("LFMM_estimated_GIF_values_", outputname, ".txt")

for (i in 1:length(pv)) {
    
    varname <- names(pv[i])
    gif <- pv[[i]]$gif
    
    cat("Estimated genomic inflation factor (gif) for:", varname, "\n", file=filename, sep=" ", append=TRUE)
    cat(gif, "\n", "--------------------------------------------------------", "\n", file=filename, sep=" ", append=TRUE)
}

# iterate over the range of GIF values and manually readjust the p-values
mygifs <- seq.int(mod.gif[1], mod.gif[2], by=by.gif)

adj.pv <- list ()
for (i in 1:length(pv)) {
    
    varname <- names(pv[i])
    
    for (j in mygifs) {
        
        gifname <- paste0("GIF_", format(j, nsmall=2))
        adj.pv[[varname]][[gifname]] <- pchisq((pv[[i]]$score[,1])^2/j, df=1, lower = FALSE)
    }
}


###### extract candidate loci ########

# create directories for result storage
lapply(prednames, function(x) if(!dir.exists(x)) dir.create(x))
for (i in prednames){
    
    if(!dir.exists(file.path(i, 'significant_SNPs'))) { dir.create(file.path(i, 'significant_SNPs')) }
    if(!dir.exists(file.path(i, 'p-value_histograms'))) { dir.create(file.path(i, 'p-value_histograms')) }
    if(!dir.exists(file.path(i, 'Manhattan_plots'))) { dir.create(file.path(i, 'Manhattan_plots')) }
    if(!dir.exists(file.path(i, 'QQ_plots'))) { dir.create(file.path(i, 'QQ_plots')) }
}

# obtain list of candidate loci by using the Benjamini-Hochberg procedure
candidates <- list()
candidate_names <- list()

for (i in 1:length(adj.pv)) {
    for (j in 1:length(adj.pv[[i]])) {
        
        varname <- names(adj.pv[i])
        gifname <- names(adj.pv[[i]][j])
        
        pvalues <- as.matrix(adj.pv[[i]][[j]])
        #pvalues[is.na(pvalues)] <- 1
        
        L = length(pvalues)
        w = which(sort(pvalues) < fdr * (1:L) / L)
        if (length(w) != 0) {
            myw <- 1:max(w)
        } else if(length(w) == 0) {
            myw <- w
        }
        
        candidates[[varname]][[gifname]] = order(pvalues)[myw]
        candidate_names[[varname]][[gifname]] <- genonames[unlist(candidates[[varname]][[gifname]])]
        
        filename = paste0("./", varname, "/significant_SNPs/", outputname, "_", varname, "_", gifname, "_LFMM_significant_loci.txt")
        write.table(candidate_names[[varname]][[gifname]], filename, sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    }
} 


###### plot results ######

# plot histogram of test significance values (it is expected to be flat, with a peak near zero)
for (i in 1:length(adj.pv)) {
    for (j in 1:length(adj.pv[[i]])) {
        
        varname <- names(adj.pv[i])
        gifname <- names(adj.pv[[i]][j])
        
        adj.pvalues <- unname(adj.pv[[i]][[j]])
        
        title_his <- paste0("./", varname, "/p-value_histograms/", outputname, "_", varname, "_", gifname, "_re-adjusted_p-value_histogram.png")
        png(filename=title_his, units="in", width=7, height=5, res=300)
        hist(adj.pvalues,
             xlab = "P-values",
             main = paste("LFMM histogram of re-adjusted p-values for", varname, gifname))
        dev.off()
        
        title_qq <- paste0("./", varname, "/QQ_plots/", outputname, "_", varname, "_", gifname, "_re-adjusted_p-value_QQplot.png")
        png(filename=title_qq, units="in", width=7, height=5, res=300)
        qqplot(rexp(length(adj.pvalues), rate = log(10)),
               -log10(adj.pvalues), 
               xlab = "Expected quantile",
               pch = 19, cex = .4)
        abline(0,1)
        title(paste("LFMM Q-Q plot of re-adjusted p-values for", varname, gifname))
        dev.off()
        
    }
} 

# plot Manhattan plot
for (i in 1:length(adj.pv)) {
    for (j in 1:length(adj.pv[[i]])) {
        
        varname <- names(adj.pv[i])
        gifname <- names(adj.pv[[i]][j])
        
        mycand <- unname(unlist(candidates[[i]][j]))
        pvalues <- unname(adj.pv[[i]][[j]])
        
        title_man <- paste0("./", varname, "/Manhattan_plots/", outputname, "_", varname, "_", gifname, "_ManhattanPlot.png")
        png(filename=title_man, units="in", width=7, height=5, res=300)
        
        plot(-log10(pvalues), 
             pch = 19, 
             cex = .5, 
             xlab = expression(bold("SNP")), ylab = expression(bold("-Log P")),
             col = ifelse(is.na(match(seq_along(-log10(pvalues)),mycand)) == TRUE,'gray50','red'))
        title(paste("LFMM Manhattan plot for", varname, gifname))
        dev.off()
    }
} 

