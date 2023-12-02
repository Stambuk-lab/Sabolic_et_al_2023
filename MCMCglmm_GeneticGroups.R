# R script for calculating heritability using Bayesian 
# MCMCglmm animal model with genetic groups.
# Based on coded example from Muff et al. (2019): Animal models with group-specific
# additive genetic variances: extending genetic group models, https://doi.org/10.1186/s12711-019-0449-7
#
# Input file should be in .txt format and look something like this
# (a set of descriptive variables followed by a set of phenotypic measurements):
#
# FAMILY    ID      MOTHER  FATHER    SEX   POP   TRAIT_1   TRAIT_2  ...  TRAIT_N
# K043K018  KK201   PK043	  PK018     M     KK    0.03950	  0.09570	 ...  0.03100
# K043K018  KM207   PM063	  PK003     F     KM   -0.00350	 -0.03210  ... -0.02270
# ...       ...     ...     ...       ...   ...   ...       ...      ...  ...
# K043K018  PK003   NA      NA        M     KK   -0.99957	 -0.64145	 ... -0.81245
# K043K018  PM063	  NA      NA        F     MM   -0.98671	 -0.64995	 ... -0.82336
# ...       ...     ...     ...       ...   ...   ...       ...      ...  ...

# needed packages:
#install.packages("MCMCglmm")
#install.packages("MasterBayes")
#install.packages("nadiv")
#install.packages("data.table")
#install.packages("pedigreemm")
library(MCMCglmm)
library(MasterBayes)
library(nadiv)
library(data.table)
library(pedigreemm)


# define input and output files:
inputfile = "siculus_common_garden_residuals.txt"       # name of the input file (the script expects .txt as an input format file)
outputfilename = "siculus_common_garden"                # prefix for the output files

# specify the pedigree data:
id = 1                                 # column containing the individual id 
mother = 3                             # column containing mothers id 
father = 2                             # column containing fathers id

# specify phenotypic traits you want to analyse:
traits = c("LwJaL",
           "LwJaO",
           "SnLgh",
           "HLgth",
           "HHght",
           "HWdth",
           "LtHip",
           "BiteF")

# specify genetic groups
gg <- c("g1"="PK",
        "g2"="PM")

# do you want to include coefficient of inbreeding as fixed effect in the model
inbreeding_coef = FALSE                 # options are TRUE or FALSE

# specify fixed and random variables to be included in the analysis:
fixedvar = c()         # specify fixed variables (if you don't want to include fixed leave EMPTY)
randomvar = c()        # specify random variables (if you don't want to include random leave EMPTY)

# define the prior distribution of residual and random effect variance (B - fixed effects variance, R - residual variance, G - random effect variances; for options see MCMCglmm package notes on CRAN)
prior <- list(R = list(V = 1, nu = 0.002),
              G = list(G1 = list(V = 1, nu = 0.002),
                       G2 = list(V = 1, nu = 0.002)))

# MCMC settings:
iteration_no = 2500000                     # define the total number of iterations 
burnin_interval = 500000                   # define the burn-in period 
thinning_interval = 100                    # define the thinning interval between consecutive observations



############################### NO NEED TO CHANGE ANYTHING UNDER THIS LINE (THEORETICALLY) #########################################


########## PREPARE THE DATASET ###########

# import data
mydata <- read.table(inputfile, header=TRUE, sep="\t", stringsAsFactors=TRUE)

# set rownames
row.names(mydata) <- mydata[[id]]

# rename the columns if necessary
colnames(mydata)[id] <- "animal"
colnames(mydata)[mother] <- "dam"
colnames(mydata)[father] <- "sire"

# extract pedigree data
mypedigreedata <- colnames(mydata[,c(id, mother, father)])
mypedigree <- mydata[, names(mydata) %in% mypedigreedata]

# arrange pedigree data so that the founders are at the top
mypedigree <- orderPed(mypedigree)

# arrange all other data according to the new pedigree order
mydata <- mydata[match(mypedigree$animal, mydata$animal),]


########## A-1 MATRIX ###########

# calculate inverse additive genetic relatedness (A-1) matrix
ainvOut <- inverseA(mypedigree)
Ainv <- ainvOut$Ainv                          #extract the A-1 matrix

# extract the coefficients of inbreeding if so specified
if (inbreeding_coef == TRUE) {
  
  f <- ainvOut$inbreeding
  mydata <- cbind(mydata, f)
}


########## GENETIC GROUPS ###########

# add genetic groups to the pedigree
my_gg <- data.frame(animal=names(gg), sire=rep(NA, length(gg)), dam=rep(NA, length(gg)))
gg_pedigree <- rbind(my_gg, mypedigree)

# replace F0's NA entries for sire and dam with genetic groups
for (i in gg) {
  
  gg_pedigree$sire[gg_pedigree$animal %like% i] <- names(gg)[gg==i]
  gg_pedigree$dam[gg_pedigree$animal %like% i] <- names(gg)[gg==i]
}

# output
write.table(gg_pedigree, paste0(outputfilename, "_GeneticGroups_pedegree.txt"), sep="\t", col.names=T, row.names=F, quote=F)

# calculate genetic groups contributions (Q) matrix
Q <- ggcontrib(gg_pedigree)

# bind Q matrix with the rest of the data
mydata <- cbind(mydata, Q)


########## PREPARE A-1 MATRIX FOR ANIMAL MODELS WITH GROUP SPECIFIC ADDITIVE VARAINCE ###########

# do the Cholesky decomposition of A-1 matrix from the pedigree
ped <- pedigree(sire=mypedigree[,"sire"], dam=mypedigree[,"dam"], label=as.character(mypedigree[,"animal"]))

# derive components T-1 and D-1 from the Cholensky decomposition
TInv <- as(ped, "sparseMatrix")
DD <- Dmat(ped)
DInv <- Diagonal(x=1/DD)

# verify that this is indeed the Cholesky decomposition by comparing the product (T-1)D-1T-1 to Ainv
invA_test<- t(TInv) %*% DInv %*% (TInv)
sum(invA_test - Ainv)

# generate A-1 for each population (group)
scaling <- list()
DD1 <- list()
DInvTilde <- list()
Ainv_gg <- list()
for (i in names(gg)) {
  
  myindex <- paste0("Index_", i)
  
  scaling[[myindex]] <- ifelse(mydata[[i]]>0,1/(mydata[[i]])^2,1e12)
  DD1[[myindex]] <- 1 - mydata[[i]]*(1-DD)
  DInvTilde[[myindex]] <- Diagonal(x=1/(DD1[[myindex]])*scaling[[myindex]])
  Ainv_gg[[myindex]] <- crossprod(sqrt(DInvTilde[[myindex]]) %*% TInv)

  # MCMCglmm requires that dimension names of the matrix have a (character) name that links columns/rows to the pedigree:
  Ainv_gg[[myindex]] <- as(Ainv_gg[[myindex]], "dgCMatrix")
  Ainv_gg[[myindex]]@Dimnames[[1L]] <- as.character(mydata$animal)
}

# prepare the dataset for calculation in MCMCglmm
for (i in names(gg)) {

  myindex <- paste0("Index_", i)
  mydata[[myindex]] <- as.character(mydata[["animal"]])
  mydata[[myindex]] <- ifelse(mydata[[i]]>0, mydata[[myindex]], NA)
}

# extract indices names
myindices <- colnames(mydata)[colnames(mydata) %like% "Index"]


########## PREPARE OUTPUT ################
filename = paste0("Bayes_MCMCglmm_results_for_", outputfilename, ".txt")

# write analysis specifications in the output file
cat("Heritability estimation using Bayesian MCMCglmm animal model", "------------------------------------------------------------------------------------", "", file=filename, sep="\n", append=TRUE)
cat("Variables: ", traits, "\n", file=filename, sep="\n", append=TRUE)
cat("Genetic groups:", gg, "\n", file=filename, sep="\n", append=TRUE)
cat("Analysis specifications:", file=filename, sep="\n", append=TRUE)
cat("Fixed variables:", fixedvar, "\n", file=filename, sep="", append=TRUE)
cat("Random variables:", randomvar, "\n", file=filename, sep=" ", append=TRUE)

cat("Prior distributions:", file=filename, sep="\n", append=TRUE)
cat(capture.output(prior), "\n", file=filename, sep="", append=TRUE)

cat("Total number of iterations:", iteration_no, "\n", file=filename, sep="", append=TRUE)
cat("Burn-in:", burnin_interval, "\n", file=filename, sep="", append=TRUE)
cat("Thinning:", thinning_interval, "\n", file=filename, sep="", append=TRUE)


########## RUN THE MODEL #################

model <- list()       # prepare a container

for (i in traits) {
      
    # set fixed and random variables
    if (is.null(fixedvar) == TRUE && inbreeding_coef == FALSE) {
      fixed <- as.formula(paste(i, "~", paste(sort(names(gg), decreasing = T), collapse = "+")))
    } else if (is.null(fixedvar) == TRUE && inbreeding_coef == TRUE) {
      fixed <- as.formula(paste(i, "~f+", paste(sort(names(gg), decreasing = T), collapse = "+")))
    } else if (is.null(fixedvar) == FALSE && inbreeding_coef == FALSE) {
      fixed <- as.formula(paste(i, "~", paste((fixedvar), collapse = "+"), "+", paste(sort(names(gg), decreasing = T), collapse = "+")))
    } else if (is.null(fixedvar) == FALSE && inbreeding_coef == TRUE) {
      fixed <- as.formula(paste(i, "~f+", paste((fixedvar), collapse = "+"), "+", paste(sort(names(gg), decreasing = T), collapse = "+")))
    } 
    
    if (is.null(randomvar) == TRUE) {
      random <- as.formula(paste("~", paste((myindices), collapse = "+")))
    } else if (is.null(randomvar) == FALSE) {
      random <- as.formula(paste("~", paste((myindices), collapse = "+"), "+", paste((randomvar), collapse = "+")))
    }
  
    # run the model
    print(i)
    k <- match(i, traits)
    
    model[[i]] <- MCMCglmm(fixed = fixed, random = random,
                           ginverse = Ainv_gg,
                           data = mydata, 
                           prior = prior,
                           singular.ok = FALSE,
                           pr = TRUE,
                           nitt = iteration_no, burn = burnin_interval, thin = thinning_interval)
}
  

########## EXPLORE THE MODEL AND ESTIMATE HERITABILITY ############

myresults <- data.frame()
herit <- list()

for (i in traits) {
  
  print(paste("Estimating heritability for:", i))
  
  # export results to output file
  cat("\n", "\n", "\n", "\n", "MCMCglmm results for: ", i, "\n", "\n",  file=filename, sep="", append=TRUE)
  cat("----------------------------------------", "\n", file=filename, sep="", append=TRUE)
  cat(capture.output(summary(model[[i]])), "\n", file=filename, sep="\n", append=TRUE)
  
  # check convergence graphically
  myplots <- paste0(outputfilename, "_trace_of_intercept_in_", i, ".pdf")
  pdf(file = myplots, width=5, height=10)
  plot(model[[i]]$Sol)
  mtext(paste(i), side = 3, line = -1.75, outer = TRUE, cex=1.25, font = 2)
  dev.off()
  
  myplots <- paste0(outputfilename, "_trace_of_variances_in_", i, ".pdf")
  pdf(file = myplots, width=10, height=5)
  plot(model[[i]]$VCV)
  mtext(paste(i), side = 3, line = -1.75, outer = TRUE, cex=1.25, font = 2)
  dev.off()
  
  # test convergence statistically (the p-values must exceed 0.05)
  cat("\n", "\n", "Heidelberger convergence diagnostic:", file=filename, sep="\n", append=TRUE)
  cat(capture.output(heidel.diag(model[[i]]$VCV)), "\n", file=filename, sep="\n", append=TRUE )
  
  # check autocorrelation
  cat("Autocorrelation in intercept:", file=filename, sep="\n", append=TRUE)
  cat(capture.output(autocorr.diag(model[[i]]$Sol)), "\n", file=filename, sep="\n", append=TRUE )
  
  cat("Autocorrelation in variances:", file=filename, sep="\n", append=TRUE)
  cat(capture.output(autocorr.diag(model[[i]]$VCV)), "\n", file=filename, sep="\n", append=TRUE )
  
  
  ########## ESTIMATE HERITABILITY ###########
  
  if (is.null(randomvar) == TRUE) {
    
    # obtain the posterior distribution of the heritability (h2 = Va/(Va+Vr))
    for (n in myindices) {
      
      herit[[i]][[n]] <- model[[i]]$VCV[, n]/(model[[i]]$VCV[, n] + model[[i]]$VCV[, "units"])
    }
    
  } else if (is.null(randomvar) == FALSE) {
    
    for (n in myindices) {
    
      # obtain the posterior distribution of the heritability (h2 = Va / (Va+Vr))
      dividend <- model[[i]]$VCV[, n]
      changing_vector <- 0
    
      for (j in 1:length(randomvar)) {
      
        changing_vector <- changing_vector + model[[i]]$VCV[, randomvar[j]]
      }
    
      divisor <- dividend + model[[i]]$VCV[, "units"] + changing_vector
    
      herit[[i]][[n]] <- dividend / divisor
    
  }
    }  
  # write heritability results in the output file
  cat("\n", "\n", "----- HERITABILITY ESTIMATION -----", "\n", file=filename, sep="", append=TRUE)
  for (n in myindices) {
    
    cat("\n", "\n", "--- ", n, " ---", "\n", file=filename, sep="", append=TRUE)
    cat("Effective size: ", effectiveSize(herit[[i]][[n]]), "\n", file=filename, sep="", append=TRUE)
    cat("Mean heritability: ", mean(herit[[i]][[n]]), "\n", file=filename, sep="", append=TRUE) 
    cat("Posterior heritability: ", posterior.mode(herit[[i]][[n]]), "\n", file=filename, sep="", append=TRUE)
    cat("95% credible interval:", HPDinterval(herit[[i]][[n]]), "\n", file=filename, sep=" ", append=TRUE) 
  }
  
  # plot
  for (n in myindices) {
  
    myplots <- paste0(outputfilename, "_heritability_distribution_for_", i, "_", n, ".png")
    png(filename = myplots, units="in", width=15, height=5, res=300)
    plot(herit[[i]][[n]])
    mtext(paste(i, "_", n), side = 3, line = -1.75, outer = TRUE, cex=1.25, font = 2)
    dev.off()
  }
  
  # export all results in a form of table
  myresults[i, "trait"] <- i
  
  for (n in myindices) {
    
    myresults[i, paste0(n, "_h2")] <- mean(herit[[i]][[n]])
    myresults[i, paste0(n, "_h2Lower")] <- HPDinterval(herit[[i]][[n]])[1]
    myresults[i, paste0(n, "_h2Upper")] <- HPDinterval(herit[[i]][[n]])[2]
    myresults[i, paste0(n, "_Va")] <- mean(model[[i]]$VCV[, n])
    myresults[i, paste0(n, "_VaLower")] <- HPDinterval(model[[i]]$VCV[, n])[1]
    myresults[i, paste0(n, "_VaUpper")] <- HPDinterval(model[[i]]$VCV[, n])[2]
    myresults[i, paste0(n, "_Vp")] <- mean(model[[i]]$VCV[, n] + model[[i]]$VCV[, "units"])
    myresults[i, paste0(n, "_VpLower")] <- HPDinterval(model[[i]]$VCV[, n] + model[[i]]$VCV[, "units"])[1]
    myresults[i, paste0(n, "_VpUpper")] <- HPDinterval(model[[i]]$VCV[, n] + model[[i]]$VCV[, "units"])[2]
  }

  myresults[i, "Ve"] <- mean(model[[i]]$VCV[, "units"])
  myresults[i, "VeLower"] <- HPDinterval(model[[i]]$VCV[, "units"])[1]
  myresults[i, "VeUpper"] <- HPDinterval(model[[i]]$VCV[, "units"])[2]
  
  if (is.null(randomvar) == FALSE) {
    
    for (j in 1:length(randomvar)) {
      
      myresults[i, randomvar[j]] <- mean(model[[i]]$VCV[, randomvar[j]])
      myresults[i, paste0(randomvar[j], "Lower")] <- HPDinterval(model[[i]]$VCV[, randomvar[j]])[1]
      myresults[i, paste0(randomvar[j], "Upper")] <- HPDinterval(model[[i]]$VCV[, randomvar[j]])[2]
    }
  }
  
  myresults[i, "DIC"] <- model[[i]]$DIC
  
  
  ########## ESTIMATE GENETIC GROUP EFFECT ###########
  
  gg_effect_col <- model[[i]][["Fixed"]][[2]]                   #calculate the number of variables to extract
  my_gg_effect <- with(model[[i]], cbind(postMode = posterior.mode(Sol[, 1:gg_effect_col]),
                                           HPDinterval(Sol[, 1:gg_effect_col])))
    
  # write results in the output file
  cat("\n", "\n", "----- GENETIC GROUPS ESTIMATION -----", "\n", file=filename, sep="", append=TRUE)
  cat("Expected difference between the mean breeding values of genetic groups: ", "\n", file=filename, sep="", append=TRUE)
  capture.output(my_gg_effect, file=filename, append=TRUE)
  
}

# export main results in a form of table
table_filename = paste0("Summarized_MCMCglmm_results_for_", outputfilename, ".txt")
write.table(myresults, table_filename, sep="\t", col.names=T, row.names=F, quote=F)

print("Analysis complete.")

