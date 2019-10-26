library(readxl)
library(dplyr)
library(PQLseq)
library(psych)
library(EMMREML)

#     Set working directory to load files -----------
#setwd()



#     Load data ---------------

# Prep cognitive and demographic data
dog_scores <- as.data.frame(read_xlsx("Supplemental_Material_dog_cognition_aging.xlsx", sheet = 9, col_names = T, range = cell_rows(2:12))) %>% 
  rename(sex = `sex*`) %>% 
  rename(reproductive.alteration = `reproductive.alteration†`)

# Prep identity by state matrix
IBS <- as.matrix(read_xlsx("Supplemental_Material_dog_cognition_aging.xlsx", sheet = 10, col_names = T, range = cell_limits(c(2,2), c(11,10))))
row.names(IBS) <-colnames(IBS)
colnames(IBS) <- NULL

# Prep Zmatrix
Zmatrix <- as.matrix(read_xlsx("Supplemental_Material_dog_cognition_aging.xlsx", sheet = 11, col_names = F, range = cell_limits(c(3,1), c(12,9))))
rownames(Zmatrix) <- seq(1:nrow(Zmatrix))



#      Model binomial cognitive abilities with PQLseq package ---------------

# Make Genetic relatedness matrix
GRM = Zmatrix%*%IBS%*%t(Zmatrix)
diag(GRM)=1 #set relatedness of the individual to 1

# Success counts
success.counts <- as.data.frame(matrix(nrow = length(dog_scores$arm.pointing), ncol = 1, data= as.numeric(dog_scores$arm.pointing)))
colnames(success.counts) <- 'success.counts'

# Total read counts
total.counts <- as.data.frame(rbind(matrix(nrow = length(dog_scores$arm.pointing), ncol = 1, data= rep(as.integer(6)))))
colnames(total.counts) <- 'total.counts'

# Covariates file
covariates <- dog_scores %>% select(sex, reproductive.alteration, mean.breed.lifespan, age) #select model covariates 

# Model matrix
mat_coeffs <- model.matrix(~ sex + reproductive.alteration + scale(mean.breed.lifespan) + poly(scale(age),2), data= covariates)

# Model binomial cognitive abilities 
out <- pqlseq(RawCountDataSet = t(success.counts), Phenotypes = mat_coeffs[,2], Covariates =  mat_coeffs[,-2], RelatednessMatrix = GRM, LibSize = t(total.counts), fit.model = 'BMM', fit.method = 'REML', fit.maxiter = 1000, fit.tol = 1e-5, filtering = FALSE, verbose = FALSE) # this line shows the set-up for using PQLseq to run a binimial mixed model with the relatedness matrix, but model with toy data will not converge



#      Generate files to model delay of gratification with EMMREML package ---------------
dog_delay_grat <- dog_scores %>% select(`delay.gratification.watching(s)`, `delay.gratification.eyes.closed(s)`, `delay.gratification.turn.back(s)`) # select only the delay of gratification tasks 

# Run Principal Component Analysis
dogs_pca <- prcomp(dog_delay_grat, center = T, scale = T) # run PCA
summary(dogs_pca) 
pca_vals <- as.data.frame(dogs_pca$x) #PCA scores for each observation
dog_delay_grat$PC1 <- pca_vals$PC1 # add PCA values to the dataframe

# Run model
model.delay.grat <- emmreml(y = dog_delay_grat$PC1, # values of 1st PCA
                        X = mat_coeffs,  # design matrix
                        Z = Zmatrix,  # matrix of random effects
                        K = IBS,  # relatedness matrix
                        varbetahat = T, varuhat = T, PEVuhat = T, test = T)

delay.grat.out <- as.data.frame(as.data.frame(model.delay.grat$pvalbeta)$none, ncol = 1); colnames(delay.grat.out) <- "pvals"  # store p-values as object
delay.grat.out$betas <- model.delay.grat$betahat
delay.grat.out$AIC <- 2*(ncol(mat_coeffs) -1) - 2*model.delay.grat$loglik # calculate AIC



