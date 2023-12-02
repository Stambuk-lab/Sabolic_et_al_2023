# R script for for running RDA analysis and 
# extracting and plotting the results.
# It will:
# import allele frequency file in Plink format, apply Hellinger transformation,
# do PCA analysis of genotype data if so specified,
# import file with additional data (i.e. phenotypic or environmental),
# run basic or partial analysis for the specified explanatory, response variables (and confounding factors),
# plot the main results and extract the results in separate .txt files,
# run permutation test for model, axis and marginal term significance, 
# and perform forward or backward explanatory variable selection for the full model if specified.

# The Plink allele frequency file looks like this:
# CHR          SNP     CLST   A1   A2      MAF    MAC  NCHROBS
# 1   2059722:34:+       BJ    G    A        0      0       40 
# 1   2059722:34:+       DU    G    A     0.84     42       50 
# 1   2059722:34:+       KL    G    A  0.06818      3       44 

# The matrix with additional data to give to model (i.e. environmental data) and constraining variables in case of partial analysis
# needs to look like this (no missing values allowed!):
# individual data -> grouping variables -> additional variables -> constraining variables (no missing values allowed):
# i.e.
#	pop1  group1    ... groupn   add1  add2  ...  addn  cons1 cons2 ... consn  
#	04	    M       ...   1      3.2   2.7   ...  4.2   3.5   2.8   ... 1.7
# 07      P       ...   2      5.7   5.8   ...  3.5   4.7   2.1   ... 2.9
#	02	    E       ...   2      0.8   1.3   ...  8.4   2.8   0.7   ... 4.3
# ...etc

# needed packages
library(reshape2)
library(vegan)
library(gridExtra)
library(ggplot2)
library(formattable)
library(RColorBrewer)


# set parameters:
genotype_file = "psic14x356x39k_m05R7r6h6DP420_1p0p.frq.strat"      # name of the Plink allele frequency file (frq.strat extension)
variable_file = "psic14_envPC123_dbMEM12.txt"                       # filename of the matrix with additional data to give to model (should be tab separated .txt)
outputname = "psic14_envPC123_dbMEM12"                              # set prefix for all output file names (no extensions)

doPCA = FALSE                    # do you want to preform PCA analysis of genotype data (options: TRUE or FALSE)
PCkeep = NULL                    # if ddPCE = TRUE set number of principal components to keep

groupingv = 2                   # set the number of grouping variables
varableno = 3                   # set the number of additional measured (environmenta) variables in variable_file

partialanalysis = TRUE         # do you want to run partial analysis (options: TRUE or FALSE)
constrainv = 2                 # if partialanalysis = TRUE set number of variables to be used as constraints

scalebiplot = 3                 # scale predictor and response scores for biplot

axis_significance = TRUE        # perform ANOVA test of RDA axes significance
margin_significance = TRUE      # perform ANOVA test of RDA marginal terms significance

selectvarforw = TRUE            # perform forward explanatory variable selection (options: TRUE or FALSE)
selectvarback = TRUE            # perform backward explanatory variable selection (options: TRUE or FALSE)


############## NO NEED TO SET ANYTHING UNDER THIS LINE ###############

options(scipen=999)                 #suppress scientific notation

###### import and prepare genotype dataset #####

print("Importing genotype dataset...")
allele_frequencies <- read.table(genotype_file, header=T, sep="", dec=".")

# extract relevant data
print("Transforming genotype dataset...")
allele_frequencies <- allele_frequencies[,c(2,3,6)]
colnames(allele_frequencies)=c("SNP","POP","MAF")

# transform to wide format
allele_frequencies <- dcast(allele_frequencies, POP~SNP)
write.table(allele_frequencies, paste0(outputname, "-freq-sites.txt"), quote=F, sep="\t", row.names=FALSE)

# apply a Hellinger transformation
geno <- allele_frequencies[,-1]       #remove pop column
geno <- decostand(geno, "hellinger")

# preform PCA on Hellinger transformed allele frequencies if so specified
if (doPCA == TRUE) {

  # run PCA
  print("Running PCA on allele frequency data...")
  geno.pca <- prcomp(geno, scale = T)
  summary(geno.pca)

  # choose PCs to keep based on broken-stick model screeplot
  png(filename=paste0(outputname,"_screeplot_geno.pca.png"), units="in", width=7, height=5, res=300)
  screeplot(geno.pca, npcs = 20, type = "lines", bstick = T)
  dev.off()
  
  # choose PCs to keep based on Kaiser-Guttman rule (those with eigenvalues greater the mean eigenvalue)
  ev <- geno.pca$sdev^2     #extract eigenvalues
  meanev <- mean(ev)        #extract mean eigenvalue
  
  filename_eig = paste0(outputname, "_eigenvalue_geno.pca.txt")
  cat(ev, "\n", "\n", file=filename_eig, sep="\t", append=TRUE)
  cat("Mean eigenvalue:", meanev, file=filename_eig, sep="", append=TRUE)
  
  # extract PC to keep
  geno.pca.axes <- geno.pca$x[,1:PCkeep]
  
  # prepare data for analysis
  genonames <- colnames(geno.pca.axes)
  responses <- geno.pca.axes
  
} else if (doPCA == FALSE) {
  
    genonames <- colnames(allele_frequencies[-1])
    responses <- geno
}


###### import additional data ######

print("Importing additional variable dataset...")

# import file with additional data to give to model
mydata <- read.table(variable_file, header = T, sep="\t")
str(mydata)


###### extract info about variables #######

print("Defining variables...")

# extract group info
colgroup = c(1:groupingv)
groupnames = colnames(mydata[colgroup])

# extract predictors
firstpredict = groupingv+1
lastpredict = groupingv+varableno
predictors = mydata[,firstpredict:lastpredict]
colpred = c(firstpredict:lastpredict)
prednames = colnames(mydata[colpred])
  

# set variables forpartial analysis if specified
if (partialanalysis == TRUE) {
  
  firstconstrain = groupingv+varableno+1
  lastconstrain = groupingv+varableno+constrainv
  constraints = mydata[,firstconstrain:lastconstrain]
  colconstrain = c(firstconstrain:lastconstrain)
  constranames = colnames(mydata[colconstrain])
  
  explainv = merge(predictors, constraints, by=0)
  explainv$Row.names <- as.numeric(explainv$Row.names)
  explainv <- explainv[order(explainv$Row.names, decreasing = FALSE),]
  rownames(explainv) <- NULL
  explainv <- explainv[,-1]
}


####### run RDA analysis ########

print("Running RDA analysis...")

if (partialanalysis == TRUE) {
  
  addq <- function(x) paste0("`", x, "`")
  form <- formula(paste("responses ~" , paste(addq(prednames), collapse = "+"), "+ Condition(", paste(constranames, collapse = "+"),")"))
  
  myrda <- rda(form, data = explainv, scale=TRUE)
  
} else if (partialanalysis == FALSE) {
  
  addq <- function(x) paste0("`", x, "`")
  form <- formula(paste("responses ~" , paste(addq(prednames), collapse = "+")))
  
  myrda <- rda(form, data=predictors, scale=TRUE)
}

summary(myrda)


####### extract results ########

print("Extracting results...")

# prepare the output file with the results of linear regression
filename_res = paste0("RDA_results_for_", outputname, ".txt")

cat("Redundancy analysis (RDA) on:", genotype_file, "and", variable_file, file=filename_res, sep="\n", append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

if (doPCA == TRUE) {
  cat("Genotype PCs are response variables. PCs used:", length(genonames), "\n", file=filename_res, sep="\n", append=TRUE)
} else if (doPCA == FALSE) {
  cat("Allele frequences are response variables. SNPs used:", length(genonames), "\n", file=filename_res, sep="\n", append=TRUE)
}

cat("Predictor variables:", prednames, "\n", file=filename_res, sep="\n", append=TRUE)

if (partialanalysis == TRUE) {
  
  cat("Constraining variables: ", constranames, "\n", file=filename_res, sep="\n", append=TRUE)
  cat("------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
} else if (partialanalysis == FALSE) {
  
  cat("No constraining variables specified for the analysis.", file=filename_res, sep="\n", append=TRUE)
  cat("------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
}

# extract main results
capture.output(myrda, file=filename_res, append=TRUE)
cat("------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

# extract scores
rdascores_sites <- scores(myrda, display="sites", choices=c(1:length(predictors*2)))
rdascores_species <- scores(myrda, display="species", choices=c(1:length(predictors*2)))
rdascores_bp <- scores(myrda, display="bp", choices=c(1:length(predictors*2)))
write.table(rdascores_sites, paste0(outputname, "_individual_RDA_scores.txt"), sep="\t") 
write.table(rdascores_species, paste0(outputname, "_response_RDA_scores.txt"), sep="\t") 
write.table(rdascores_bp, paste0(outputname, "_predictor_RDA_scores.txt"), sep="\t") 

# extract regression coefficients
mylrcoeff <- coef(myrda)
write.table(mylrcoeff, paste0(outputname, "_RDA_linear_regression_coefficients.txt"), sep="\t")

# variance explained
myeigenvalues <- summary(eigenvals(myrda, model = "constrained"))
myxlab <- paste0("RDA1 (", percent(myeigenvalues[2,1]), ")")
myylab <- paste0("RDA2 (", percent(myeigenvalues[2,2]), ")")
capture.output(myeigenvalues, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

# extract R2 values
myR2 <- RsquareAdj(myrda)$r.squared
myR2adj <- RsquareAdj(myrda)$adj.r.squared
cat("R2:", myR2, "\n", file=filename_res, sep=" ", append=TRUE)
cat("R2 adjusted:", myR2adj, "\n", file=filename_res, sep=" ", append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

# permutation tests for RDA results
print("Running permutation tests...")

print("Checking the significance of the model...")
myaov <- anova.cca(myrda, parallel = getOption("mc.cores", 2))

print("Extracting results...")
cat("ANOVA permutation test for RDA model:", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(myaov, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)

png("AOV_model.png")
a1 <- tableGrob(myaov)
grid.arrange(a1)
dev.off()

if (axis_significance == TRUE) {
  
  print("Checking the significance of RDA axes...")
  myaovax <- anova.cca(myrda, by='axis', parallel = getOption("mc.cores", 2))
  
  print("Extracting results...")
  cat("ANOVA significance for each constrained axis:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(myaovax, file=filename_res, append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
  png("AOV_axis.png")
  a2 <- tableGrob(myaovax)
  grid.arrange(a2)
  dev.off()
}

if (margin_significance == TRUE) {
  
  print("Checking the significance predictors...")
  myaovmg <- anova.cca(myrda, by='margin', parallel = getOption("mc.cores", 2))
  
  print("Extracting results...")
  cat("ANOVA significance for marginal effect of predictors:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(myaovmg, file=filename_res, append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
  png("AOV_marginal_effect.png")
  a3 <- tableGrob(myaovmg)
  grid.arrange(a3)
  dev.off()
}

# goodness of fit
mygoodnessoffit <- goodness(myrda)
write.table(mygoodnessoffit, paste0(outputname, "_RDA_goodness_of_fit.txt"), sep="\t")

# variance inflation factors
myvif <- vif.cca(myrda)
  
cat("Variance inflation factors (VIF):", "\n", file=filename_res, sep="\n", append=TRUE)
capture.output(myvif, file=filename_res, append=TRUE)
cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)


##### explanitory variable selection ######

# forward explanatory variable selection
if (selectvarforw == TRUE) {
  
  print("Running forward explanatory variable selection...")
  
  if (partialanalysis == TRUE) {
    
    
    mod0 <- rda(responses~1, data = explainv)
    
  } else if (partialanalysis == FALSE) {
    
    mod0 <- rda(responses~1, data = predictors)
  }
  
  step.forward <- ordistep(mod0, scope=formula(myrda), direction="forward", permutations=how(nperm=999))
  myR2adj_ordi <- (RsquareAdj(step.forward))$adj.r.squared
  
  cat("Forward ordistep selection of explanatory variables:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(step.forward, file=filename_res, append=TRUE)
  cat("\n", "\n", file=filename_res, sep=" ", append=TRUE)
  cat("Forward ordistep selection R2 adjusted:", myR2adj_ordi, file=filename_res, sep=" ", append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
  step2.forward <- ordiR2step(mod0, scope=formula(myrda), direction="forward", R2scope=TRUE, permutations=how(nperm=999))
  myR2adj_ordi <- (RsquareAdj(step2.forward))$adj.r.squared
  
  cat("Forward ordiR2step selection of explanatory variables:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(step.forward, file=filename_res, append=TRUE)
  cat("\n", "\n", file=filename_res, sep=" ", append=TRUE)
  cat("Forward ordiR2step selection R2 adjusted:", myR2adj_ordi, file=filename_res, sep=" ", append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
}


# backward explanatory variable selection
if (selectvarback == TRUE) {
  
  print("Running backward explanatory variable selection...")
  
  myordi <- ordistep(myrda, permutations=how(nperm=999))
  myR2adj_ordi<- (RsquareAdj(myordi))$adj.r.squared
  
  cat("Backward ordistep selection of explanatory variables:", "\n", file=filename_res, sep="\n", append=TRUE)
  capture.output(myordi, file=filename_res, append=TRUE)
  cat("\n", "\n", file=filename_res, sep=" ", append=TRUE)
  cat("Backward selection R2 adjusted:", myR2adj_ordi, file=filename_res, sep=" ", append=TRUE)
  cat("\n", "------------------------------------------------------------------------------------", "\n", file=filename_res, sep="\n", append=TRUE)
  
}


##### plot RDA ######

print("Plotting RDA...")

# plot R2 scores
myrdaR2 <- data.frame(score=c("R2", "R2adj"), value=c(myR2, myR2adj))

title_rplot <- paste0(outputname, "_R2.png")
png(filename=title_rplot, units="in", width=5, height=5, res=300)
gp <- ggplot(data=myrdaR2, aes(x=score, y=value)) +
  geom_bar(stat = "identity", color="black", fill=c("forestgreen","darkorange")) +
  scale_x_discrete(limits=myrdaR2$score, labels = c("R2","R2adj")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(0,1)) +
  theme_bw() +
  labs(title=expression(paste(R^2," and ",R^2," adjusted scores"))) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(size=12, face="bold"),
        plot.title = element_text(size=18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(gp)

dev.off()


# plot rest
if (length(myeigenvalues)>3) {
  
  # plot eigenvalues
  title_eig <- paste0(outputname,"_RDA_eigenvalues", ".png")
  png(filename=title_eig, units="in", width=5, height=5, res=300)
  screeplot(myrda, main="RDA eigenvalues")
  dev.off()
  
  # plot triplot
  title_scattplot <- paste0(outputname,"_RDA_triplot", ".png")
  png(filename=title_scattplot, units="in", width=7, height=7, res=300)
  plot(myrda, scaling=scalebiplot, type="n", xlab=myxlab, ylab=myylab, font.lab=2)
  points(myrda, scaling=scalebiplot, display="sites", pch=20, cex=1, col="gray55")                 #add individuals scores
  points(myrda, scaling=scalebiplot, display="species", pch=20, cex=1.5, col="red")                #add phenotypic measures 
  text(myrda, scaling=scalebiplot, display="species", col="red", cex=0.8, font=2, pos=3)           #add phenotypic measures
  text(myrda, scaling=scalebiplot, display="bp", col="black", cex=1, font=2)                       #add predictor biplot scores
  dev.off()
  
  # plot grouping scores 
  for ( i in groupnames) {
    
    # extract group levels
    mydata[,i] <- as.factor(mydata[,i])
    tags <- levels(mydata[[i]])
    numgroups <- length(tags)
    
    # set colours
    if (numgroups<=4)	{
      mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3")
      myshapes=c(22,21,23,24)
    } else if (numgroups>4 & numgroups<=11) {	
      mypalette <- c("darkorange3", "darkgreen", "orange", "chartreuse3", "black", "red1", "goldenrod1", "purple3", "blue", "magenta", "cyan")
      myshapes=c(21:25,21:25,21)
    } else if (numgroups==12) {	
      mypalette <- brewer.pal(12,"Paired")
      myshapes=c(21:25,21:25,21:22)
    } else if (numgroups==14) {	
      mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(14)
      myshapes=c(21:25,21:25,21:24)
    } else if (numgroups==20) {	
      mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(20)
      myshapes=c(21:25,21:25,21:25,21:25)
    } else if (numgroups==26) {	
      mypalette <- colorRampPalette(brewer.pal(name="Paired", n = 12))(26)
      myshapes=c(21:25,21:25,21:25,21:25,21:25,21)
    }
    
    # plot
    title_scattplot <- paste0(outputname, "_", i, "_RDA", ".png")
    png(filename=title_scattplot, units="in", width=7, height=5, res=300)
    plot(myrda, type = "n", xlab=myxlab, ylab=myylab, font.lab=2)
    with(mydata, points(myrda, display = "sites", col = mypalette[mydata[[i]]], pch = 21, bg=mypalette[mydata[[i]]]))
    with(mydata, legend("topleft", legend = levels(mydata[[i]]), bty = "n", col = mypalette, pch = 21, pt.bg = mypalette, text.font=2, cex=0.55))
    
    dev.off()
    
  }
  
} else if (length(myeigenvalues)<=3) {
  
  print("Can't plot, only one RDA eigenvector obtained.")
}


