# R script to run PCA analysis on genomic datasets.
# The input file needs to be in VCF format, 
# and a corresponding popmap file should be provided.

# needed packages
#install.packages("vcfR")
#install.packages("adegenet")
#install.packages("ggfortify")
library(vcfR)
library(adegenet)
library(ggfortify)


# set needed parameters

# set input VCF file name
myvcf = "sicula585x39k_m05R7r6h6DP420_1p0p.vcf"

# set corresponding popmap file name
mypopmap = "popmap_sicula585x39k_m05R7r6h6DP420_1p0p"

# set the number of PCs to keep
keeppcs = 5

# set output file name suffix
output_name = "sicula585x39k_m05R7r6h6DP420_1p0p"



#############################################################################################################


# import VCF file
rawvcf <- read.vcfR(myvcf)

# transform to GENIND format
rawgenind <- vcfR2genind(rawvcf)

# import population information
pop <- read.table(mypopmap, header=F, sep="\t", stringsAsFactors = TRUE)

# add to GENIND file
rawgenind@pop <- pop[,2]

# scale GENIND data
PCAmatrix <- scaleGen(rawgenind)   

# run PCA
PCAraw <- prcomp(PCAmatrix)

# extract explained variance
prop.pca = PCAraw$sdev^2/sum(PCAraw$sdev^2)
mypc_labels <- list()
for (i in 1:keeppcs) {
  mypc <- paste0("PC", i)
  mypc_labels[[mypc]] <- paste0(mypc, " (", round(100*(prop.pca[i]), digits=2), "%)")
}

# transform matrix to dataframe and add population information
PCA_matrix_data <- as.data.frame(PCAmatrix)
PCA_matrix_data$Pops <- as.factor(pop[,2])

# set colors
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

# plot PCA scatter plots
for (i in 1:keeppcs) {
  
  first_pc = i
  second_pc = i+1
  
  if (second_pc <= keeppcs) {
    
    # plot PCA
    plotname = paste0("PCA_plot_PC", first_pc, "-PC", second_pc, "_", output_name, ".png")
    png(filename=plotname, units="in", width=8, height=5, res=600)
    print(
      autoplot(PCAraw, x=first_pc, y=second_pc, data=PCA_matrix_data) +
        geom_point(aes(fill=as.factor(PCA_matrix_data$Pops)), color="gray30", shape=21, size=3) +
        scale_fill_manual(values=mypallete) +
        labs(x = mypc_labels[[paste0("PC", first_pc)]], y = mypc_labels[[paste0("PC", second_pc)]]) +
        theme_bw() +
        theme(legend.title=element_blank(), legend.text=element_text(size=10, face="bold"),
              axis.title=element_text(size=14, face="bold"), axis.text=element_text(size=10))
    )
    dev.off()
  }
}

# plot PCA scatter plot with individual information
indplotname = paste0("PCA_individuals_PC1-PC2_", output_name, ".png")
png(filename=indplotname, units="in", width=8, height=5, res=600)
ggplot(PCAraw$x, aes(x=PC1, y=PC2)) +
  geom_point(size = 0.5) + 
  geom_text(aes(label=row.names(PCAraw$x)), size=1.5, hjust=0.5, vjust=0, nudge_y=0.5) +
  labs(x = mypc_labels[["PC1"]], y = mypc_labels[["PC2"]]) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.text=element_text(size=10, face="bold"),
        axis.title=element_text(size=14, face="bold"), axis.text=element_text(size=10))
dev.off()

# plot eigenvalues
prop.pca2 <- prop.pca[1:10]
names(prop.pca2) <- 1:10
png(filename=paste0("explained_variance", output_name, ".png"), units="in", width=5, height=5, res=300)
barplot(prop.pca2, xlab="Principal component", ylab="Proportion of Explained Variance", cex.lab=1, cex.names=0.8, ylim=c(0,0.04))
dev.off()

# extract results summary
out <- capture.output(summary(PCAraw))	
cat(out, file=paste0("PCAsummary_", output_name, ".txt"), sep="\n", quote=FALSE)

# exctract loadings
loadings <- PCAraw$rotation
write.table(format(loadings, digits=4),file=paste0("PCAloadings_", output_name, ".txt"), sep="\t", quote=FALSE)

# ectract scores
pcascores <- PCAraw$x
write.table(pcascores,file=paste0("PCAscores_", output_name, ".txt"), sep="\t", quote=FALSE)

