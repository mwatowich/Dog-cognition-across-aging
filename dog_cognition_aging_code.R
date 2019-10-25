library(readxl)
library(dplyr)
library(PQLseq)
library(psych)
library(EMMREML)

# Set working directory -----------
#setwd()

#     Load data ---------------
# Prep Dognition cognitive and demographic data
dog_scores <- as.data.frame(read_xlsx("Dognition_Supplemental_Material.xlsx", sheet = 9, col_names = F))
dog_scores <- dog_scores[-1,]
colnames(dog_scores) <- dog_scores[1,]
dog_scores <- dog_scores[-1,]

# Prep identity by state matrix
dog_names <- c("Shih_Tzu", "Labrador_Retriever", "Golden_Retriever", "Dalmatian", "Shetland_Sheepdog", "Siberian_Husky", "Vizsla", "Pug", "Boxer")
ibs_tmp <- as.matrix(read_xlsx("Dognition_Supplemental_Material.xlsx", sheet = 10, col_names = F))
ibs <- ibs_tmp[-c(1:2),-1]
IBS <- matrix(as.numeric(ibs[1:nrow(ibs), 1:ncol(ibs)]), ncol = ncol(ibs), nrow = nrow(ibs))
rownames(IBS) <- dog_names

# Prep Zmatrix
Zmatrix_tmp <- as.matrix(read_xlsx("Dognition_Supplemental_Material.xlsx", sheet = 11, col_names = F))
Zmatrix <- Zmatrix_tmp[-c(1:2),]
Zmatrix <- apply(Zmatrix, 2, as.numeric)
rownames(Zmatrix) <- seq(1:nrow(Zmatrix))


#      Generate files to model binomial cognitive abilities with PQLseq package ---------------

# Make Genetic relatedness matrix
GRM = Zmatrix%*%IBS%*%t(Zmatrix)
diag(GRM)=1 #set relatedness of the individual to 1

# Success counts
success.counts <- as.data.frame(matrix(nrow = length(dog_scores$arm.pointing), ncol = 1, data= as.integer(dog_scores$arm.pointing)))
colnames(success.counts) <- 'success.counts'

# Total read counts
total.counts <- as.data.frame(rbind(matrix(nrow = length(dog_scores$arm.pointing), ncol = 1, data= rep(as.integer(6)))))
colnames(total.counts) <- 'total.counts'

# Covariates file
dog_scores$sex <- ifelse(dog_scores$sex == "male", 1, 0)
covariates <- dog_scores %>% select(sex, reproductive.alteration, mean.breed.lifespan, age) #columns we want in pqlseq
covariates <- as.data.frame(apply(covariates, 2, as.numeric))
colnames(covariates) <- c('sex', 'reproductive.alteration', 'mean.breed.lifespan', 'age') 
covariates$scale.age <- as.numeric(scale(covariates$age)) #scale variables 

# Model matrix
mat_coeffs <- model.matrix(~ sex + reproductive.alteration + scale(mean.breed.lifespan) + poly(scale.age,2), data= covariates)


#      Model binomial cognitive abilities with PQLseq package ---------------
out <- pqlseq(RawCountDataSet = t(success.counts), Phenotypes = mat_coeffs[,2], Covariates =  mat_coeffs[,-2], RelatednessMatrix = GRM, LibSize = t(total.counts), fit.model = 'BMM', fit.method = 'REML', fit.maxiter = 1000, fit.tol = 1e-5, filtering = FALSE, verbose = FALSE) # this line shows the set-up for using PQLseq to run a binimial mixed model with the relatedness matrix, but model with toy data will not converge

out$AIC <- 2*(ncol(mat_coeffs) -1) - 2*out$loglik #calculate AIC


#      Generate files to model delay of gratification with EMMREML package ---------------
dog_delay_grat <- dog_scores %>% select(delay.gratification.watching, delay.gratification.eyes.closed, delay.gratification.turn.back) # select only the delay of gratification tasks 
dog_delay_grat <- as.data.frame(apply(dog_delay_grat, 2, as.numeric)) # ensure these are numeric values

# Run Principal Component Analysis
dogs_pca <- prcomp(dog_delay_grat, center = T, scale = T) # run PCA
summary(dogs_pca) #variance explained
fviz_eig(dogs_pca) #scree plot
head(dogs_pca$rotation) #PCA loadings 
pca_vals <- as.data.frame(dogs_pca$x) #PCA scores for each observation
dog_delay_grat$PC1 <- pca_vals$PC1 # add PCA values to the dataframe. These are the outcome variable in the model below

# Run model
model.delay.grat <- emmreml(y = dog_delay_grat$PC1, # values of 1st PCA
                        X = mat_coeffs,  # design matrix
                        Z = Zmatrix,  # matrix of random effects
                        K = IBS,  # relatedness matrix
                        varbetahat = T, varuhat = T, PEVuhat = T, test = T)

delay.grat.out <- as.data.frame(as.data.frame(model.delay.grat$pvalbeta)$none, ncol = 1); colnames(delay.grat.out) <- "pvals"  # store p-values as object
delay.grat.out$betas <- model.delay.grat$betahat
delay.grat.out$AIC <- 2*(ncol(mat_coeffs) -1) - 2*model.delay.grat$loglik # calculate AIC



