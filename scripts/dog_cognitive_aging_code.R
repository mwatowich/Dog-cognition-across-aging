#!/usr/bin/env Rscript

#     Load libraries
library(tidyverse)
library(PQLseq)
library(psych)
library(EMMREML)


#     Load data ---------------
# Download the toy IBD matrix, Z matrix, and data from the Github page mwatowich/Dog-cognition-across-aging/scripts.

# Load cognitive and demographic data
dog_scores <- read_csv(file = "example_data_dog_cognitive_aging.csv", 
                       col_names = T, n_max = 10) %>% 
  rename(sex = `sex*`) %>% 
  rename(reproductive.alteration = `reproductive.alteration†`) %>% 
  as.data.frame()

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


#      Model binomial cognitive tasks with PQLseq package ---------------
# Genetic relatedness matrix
GRM = Zmatrix%*%IBD%*%t(Zmatrix) # multiply Zmatrix and IBD matrices to create genetic relatedness matrix
diag(GRM)=1 # set relatedness of the individual to 1

# Success counts
success.counts <- as.data.frame(matrix(nrow = length(dog_scores$arm.pointing), 
                                       ncol = 1, 
                                       data= as.numeric(dog_scores$arm.pointing)))
colnames(success.counts) <- c('success.counts') # change column names 

# Total read counts
total.counts <- as.data.frame(rbind(matrix(nrow = length(dog_scores$arm.pointing), 
                                           ncol = 1, 
                                           data= rep(as.integer(6)))))
colnames(total.counts) <- c('total.counts') # change column names 

# Covariates file
covariates <- dog_scores %>% 
  select(sex, 
         reproductive.alteration, 
         mean.breed.lifespan, 
         age) #select model covariates 

# Model matrix
mat_coeffs <- model.matrix(~ sex + 
                             reproductive.alteration + 
                             scale(mean.breed.lifespan) + 
                             poly(scale(age),2), 
                           data= covariates)

# Model binomial cognitive performance 
# This shows the set-up for using PQLseq to run a binimial mixed model with the relatedness matrix 
# Please note that model with toy data will not converge 
out <- pqlseq(RawCountDataSet = t(success.counts), 
              Phenotypes = mat_coeffs[,2], 
              Covariates =  mat_coeffs[,-2], 
              RelatednessMatrix = GRM, 
              LibSize = t(total.counts), 
              fit.model = 'BMM', 
              fit.method = 'REML') 


#      Model delay of gratification (or eye contact) with EMMREML ---------------
# Select only the delay of gratification tasks 
dog_delay_grat <- dog_scores %>% 
select(`delay.gratification.watching(s)`, 
       `delay.gratification.eyes.closed(s)`, 
       `delay.gratification.turn.back(s)`) 

# Principal Component Analysis
dogs_pca <- prcomp(dog_delay_grat, center = T, scale = T) # run PCA
pca_vals <- as.data.frame(dogs_pca$x) #PCA scores for each observation
dog_delay_grat$PC1 <- pca_vals$PC1 # add PCA values to the dataframe

# Run model 
model.delay.grat <- emmreml(y = dog_delay_grat$PC1, # values of 1st PCA
                        X = mat_coeffs,  # design matrix
                        Z = Zmatrix,  # matrix of random effects
                        K = IBD,  # relatedness matrix
                        varbetahat = T, varuhat = T, PEVuhat = T, test = T)

delay.grat.out <- cbind(model.delay.grat$betahat, model.delay.grat$pvalbeta[, "none"])
colnames(delay.grat.out) <- c("p-values", "betas")
