# Rscript to run PCAdapt genome scan to detect loci under selection.
# Input file is in VCF format.
# Popmap (two columns with individual and population IDs),
# and locilist files are required (no headers allowed in either file).

# needed packages:
#install.packages("pcadapt")
#install.packages("ggplot2")
library(pcadapt)
library(ggplot2)


# set input file names:
vcf_filename = "NPK-NPM82x10k_m05R7r6h6DP420_1p0p.vcf"
popmap_filename="popmap_NPK-NPM82x10k_m05R7r6h6DP420_1p0p"
locilist_filename="locilist_NPK-NPM82x10k_m05R7r6h6DP420_1p0p"

# set the prefix for the output files:
outputname = "NPK-NPM82x10k_m05R7r6h6DP420_1p0p_K2_T0.05"

# set the maximal number of K values you want to test
max_K = 10

# set the chosen K value you want to use to detect outliers
choosen_K = 1

# set the minimum allele frequency threshold
min_maf_threshold = 0.05

# set outlier p-value significance threshold
sig_threshold = 0.05


################################################################################################################

# import VCF file
vcf_file <- read.pcadapt(vcf_filename, type="vcf")

# import popmap file and edit colnames
pop_map <- read.table(popmap_filename, header=FALSE, sep="\t")
colnames(pop_map) <- c("IND", "POP")

# import locilist file and edit colnames
locilist <- read.table(locilist_filename, header=FALSE, sep="\t")
colnames(locilist) <- "snp"

# run a trial PCAdapt analysis with large number of K values
pcadapt_trial <- pcadapt(vcf_file, K=max_K, min.maf=min_maf_threshold) 

# plot screeplot
png(filename=paste0(outputname, "PCAdapt_screeplot.png"), units="in", width=8, height=5, res=600)
print(
  plot(pcadapt_trial, option = "screeplot", plt.pkg = "ggplot"))
dev.off()

# set colors for scatterplots
mypallete <- c("dodgerblue2",
               "green4",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black",
               "skyblue2",
               "#FB9A99", # lt pink
               "palegreen2",
               "#E31A1C", # red
               "#CAB2D6", # lt purple
               "#FDBF6F", # lt orange
               "gray40",
               "khaki2",
               "maroon",
               "orchid1",
               "deeppink1",
               "blue1",
               "steelblue4",
               "gold1")

# plot scatter plot for PC1 and PC2
png(filename=paste0(outputname, "PCAdapt_scatteplot_PC1-PC2.png"), units="in", width=8, height=5, res=600)
print(
  plot(pcadapt_trial, option = "scores", pop=pop_map$POP, i=1, j=2, col=mypallete, plt.pkg="ggplot"))
dev.off()

# plot scatter plot for  PC3 and PC4
png(filename=paste0(outputname, "PCAdapt_scatteplot_PC3-PC4.png"), units="in", width=8, height=5, res=600)
print(
  plot(pcadapt_trial, option = "scores", pop=pop_map$POP, i=3, j=4, col=mypallete, plt.pkg = "ggplot"))
dev.off()


# run PCAdapt again with the choosen K value
pcadapt_results <- pcadapt(vcf_file, K=choosen_K, min.maf=min_maf_threshold) 

# plot Manhattan plot
png(filename=paste0(outputname, "PCAdapt_Manhattan.png"), units="in", width=8, height=5, res=600)
print(
  plot(pcadapt_results, option = "manhattan"))
dev.off()

# plot Q-Q plot
png(filename=paste0(outputname, "PCAdapt_QQplot.png"), units="in", width=8, height=5, res=600)
print(
  plot(pcadapt_results, option = "qqplot"))
dev.off()

# plot histograms of the p-values
png(filename=paste0(outputname, "PCAdapt_p-value_histogram.png"), units="in", width=8, height=5, res=600)
print(
  hist(pcadapt_results$pvalues, xlab="p-values", main=NULL, breaks=50, col="gray60"))
dev.off()

# extract unadjusted PCAdapt significant loci
significant_snps <- locilist[which(pcadapt_results$pvalues < sig_threshold),] 

# export unadjusted PCAdapt significant loci
cat(significant_snps, file=paste0(outputname, "PCAdapt_unadjusted_outliers.txt"), sep="\n")

# extract uadjusted p-values into separate dataframe
unadjusted_pvalues <- as.data.frame(cbind(c(1:length(locilist$snp)), -log10(pcadapt_results$pvalues)))
names(unadjusted_pvalues) <- c("Loci", "pvalue")

# plot Manhattan plot with significant loci marked in red
png(filename=paste0(outputname, "PCAdapt_Manhattan_unadjusted_pvalues.png"), units="in", width=8, height=5, res=600)
print(
  ggplot(unadjusted_pvalues, aes(x=Loci, y=pvalue, col=ifelse(pvalue < -log10(sig_threshold),'gray50','red'))) + 
    geom_point(cex=2) +
    scale_colour_identity() +
    theme_bw() + theme(legend.position="none", 
                       axis.title.x=element_text(size=14,face="bold"), axis.title.y=element_text(size=14, vjust=1),
                       axis.text=element_text(size=10)) +
    geom_line(aes(y=-log10(sig_threshold), colour="black"), linetype="dashed", size=0.8) +
    labs(y=expression(bold("-log10(p-value)"))))
dev.off()

# run Benjamini-Hochberg procedure to adjust p-values 
pcadapt_results_padj <- p.adjust(pcadapt_results$pvalues,method="BH")

# extract BH-adjusted PCAdapt significant loci
significant_snps_padj <- locilist[which(pcadapt_results_padj < sig_threshold),]

# export BH-adjusted PCAdapt significant loci
cat(significant_snps_padj, file=paste0(outputname, "PCAdapt_BH-adjusted_outliers.txt"), sep="\n")

# extract BH-adjusted p-values into separate dataframe
adjusted_pvalues <- as.data.frame(cbind(c(1:length(locilist$snp)), -log10(pcadapt_results_padj)))
names(adjusted_pvalues) <- c("Loci", "pvalue")

# plot Manhattan plot with significant loci marked in red
png(filename=paste0(outputname, "PCAdapt_Manhattan_BH-adjusted_pvalues.png"), units="in", width=8, height=5, res=600)
print(
  ggplot(adjusted_pvalues, aes(x=Loci, y=pvalue, col=ifelse(pvalue < -log10(sig_threshold),'gray50','red'))) + 
    geom_point(cex=2) +
    scale_colour_identity() +
    theme_bw() + theme(legend.position="none", 
                       axis.title.x=element_text(size=14,face="bold"), axis.title.y=element_text(size=14, vjust=1),
                       axis.text=element_text(size=10)) +
    geom_line(aes(y=-log10(sig_threshold), colour="black"), linetype="dashed", size=0.8) +
    labs(y=expression(bold("-log10(p-value)"))))
dev.off()



