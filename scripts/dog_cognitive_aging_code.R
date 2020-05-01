#!/usr/bin/env Rscript

#     Load libraries
library(tidyverse)
library(PQLseq)
library(psych)
library(EMMREML)


#     Load data ---------------
# Download the example IBD matrix, Z matrix, and data from the Github page mwatowich/Dog-cognition-across-aging 

# Load cognitive and demographic data
dog_scores <- read_csv(file = "example_data_dog_cognitive_aging.csv", 
                       col_names = T, n_max = 10) %>% as.data.frame()
colnames(dog_scores)[c(4,5)] <- c("sex", "reproductive.alteration")

# Load identity by decent matrix
IBD <- read_delim(file = "example_IBDmatrix_dog_cognitive_aging.txt", 
                  delim = "\t", 
                  col_names = T) %>% 
  column_to_rownames("X1") %>% 
  as.matrix()

# Load Zmatrix
Zmatrix <- as.matrix(read_csv("example_Zmatrix_dog_cognitive_aging.csv", 
                              col_names = T))
rownames(Zmatrix) <- seq(1:nrow(Zmatrix))


#      Make model matrices to test equations 1-4 from Watowich et al. 2020 ---------------
# Model matrix with all covariates and age as a linear function (Eq 1)
linear_cov <- model.matrix(~ sex + 
                             reproductive.alteration +
                             scale(mean.breed.lifespan) + 
                             scale(age),
                           data= dog_scores)

# Model matrix with age as a logarithmic function (Eq 2)
log_cov <- model.matrix(~ sex + 
                          reproductive.alteration +
                          scale(mean.breed.lifespan) + 
                          log(2 + scale(age)), 
                        data= dog_scores)

# Model matrix of the truncation hypothesis (additive; Eq 3)
truncation_cov <- model.matrix(~ sex + 
                                 reproductive.alteration +
                                 scale(mean.breed.lifespan) + 
                                 poly(scale(age),2),
                               data= dog_scores)

# Model matrix of the compression hypothesis (interactive; Eq 4)
compression_cov <- model.matrix(~ sex + 
                                  reproductive.alteration +
                                  scale(mean.breed.lifespan)*poly(scale(age),2), 
                                data= dog_scores)


#      Model binomial cognitive tasks with PQLseq package ---------------
# Genetic relatedness matrix
GRM = Zmatrix%*%IBD%*%t(Zmatrix)                # multiply Zmatrix and IBD matrices to create genetic relatedness matrix
diag(GRM)=1                                     # set relatedness of the individual to 1

# Success counts
success.counts <- as.data.frame(matrix(nrow = length(dog_scores$arm.pointing), 
                                       ncol = 1, 
                                       data= as.numeric(dog_scores$arm.pointing)))
colnames(success.counts) <- c('success.counts') # change column names 

# Total read counts
total.counts <- as.data.frame(rbind(matrix(nrow = length(dog_scores$arm.pointing), 
                                           ncol = 1, 
                                           data= rep(as.integer(6)))))
colnames(total.counts) <- c('total.counts')     # change column names 


# Model binomial tasks of cognitive performance
# This shows the set-up for using PQLseq to run a binimial mixed model with a genetic relatedness matrix 
# Please note that models with example data will not converge 

# Example: Model of the truncation hypothesis 
truncation_out <- pqlseq(RawCountDataSet = t(success.counts), 
                         Phenotypes = truncation_cov[,2],    # choose one covariate to set as the "phenotype" 
                         Covariates =  truncation_cov[,-2],  # set all other covariates as "covariates" 
                         RelatednessMatrix = GRM, 
                         LibSize = t(total.counts),
                         fit.model = 'BMM', 
                         fit.method = 'REML') 

# Model the linear, log, truncation, and compression models (Eq 1-4) as shown above 

# Calculate AIC from loglik output, compare AIC of the linear, log, and quadratic models (Eq 1-3) 
# Calculate AIC from loglik output, compare AIC of the truncation and compression models (Eq 3-4)


#      Model delay of gratification (or eye contact) with EMMREML ---------------
# Select only the delay of gratification tasks 
dog_delay_grat <- dog_scores %>% 
  dplyr::select(`delay.gratification.watching(s)`, 
       `delay.gratification.eyes.closed(s)`, 
       `delay.gratification.turn.back(s)`) 

# Principal Component Analysis
dogs_pca <- prcomp(dog_delay_grat, center = T, scale = T)   # perform PCA
pca_vals <- as.data.frame(dogs_pca$x)                       # PCA scores for each observation
dog_delay_grat$PC1 <- pca_vals$PC1                          # add PC1 values to the dataframe


# Example: Model of the truncation hypothesis 
truncation_delay_grat <- emmreml(y = dog_delay_grat$PC1,    # values of 1st PCA
                        X = truncation_cov,                 # covariates 
                        Z = Zmatrix,                        # matrix of random effects
                        K = IBD,                            # relatedness matrix
                        varbetahat = T, varuhat = T, PEVuhat = T, test = T)

truncation_delay_grat <- cbind(truncation_delay_grat$betahat, 
                               truncation_delay_grat$varbetahat, 
                               truncation_delay_grat$pvalbeta[, "none"], 
                               truncation_delay_grat$loglik)

colnames(truncation_delay_grat) <- c("betas", "se", "p-values", "loglik")

# Model the linear, log, truncation, and compression models (Eq 1-4) as shown above 

# Calculate AIC from loglik output, compare AIC of the linear, log, and quadratic models (Eq 1-3) 
# Calculate AIC from loglik output, compare AIC of the truncation and compression models (Eq 3-4)
