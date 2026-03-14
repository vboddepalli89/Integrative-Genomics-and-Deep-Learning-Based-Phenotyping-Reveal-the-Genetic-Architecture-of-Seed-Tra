# Define the file path
file_path <- r"(C:\Users\naresh89\Desktop\My projects\gwas\Seed image analysis\GWAS_Weilv9\SVEN\sven_weilv9\seed_numerical.txt)"
gz_file_path <- paste0(file_path, ".gz")

# Use gzfile to compress
file_content <- readLines(file_path)
gz_con <- gzfile(gz_file_path, "w")
writeLines(file_content, gz_con)
close(gz_con)

cat("File compressed to:", gz_file_path, "\n")

library(Matrix) # To create sparse matrices
getwd()
# MLScoreTraits_imputed.txt.gz is obtained from TASSEL after
# first running a kNN imputation and then a simple average
# imputation. Saved as numerical genotype. Then compressed
# using gzip to save space on disk.

# read the numeric genotype output from TASSEL
snp = read.table("seed_numerical.txt.gz",header=T)
dim(snp)
snp$X.Marker.
# The first column contains the name of the plants. Store the SNP values in X:
X = as.matrix(snp[,-1])


# Since major alleles are coded as 1 and minor alleles are 0
# reverse the coding to have more zeros (which can be ignored)
# This also does not affect the model except for a sign change of the coefficients.
X = 1 - X

# Check: the following should be much less than 0.5 for sparse format to be useful
# If you find the answer is NA, there's probably still SNPs with missing alleles values.
print(mean(X != 0))

# Now find columns with MAF > 0.05
mx = colMeans(X) > 0.05

# Convert to sparse format, first filetering out columns wiht MAF < 0.05
Z = Matrix(X[,mx],sparse = T)

# assign the row names for later reference....
# we need these to make sure the phenotype values are correctly
# matched with the genotypes.
rownames(Z) <- snp$X.Marker.
saveRDS(Z,file="zmat_sparse.rds")

##################
# Check out SNP_to_X.r if you have not already
#install.packages("bravo")
#install.packages("doParallel")
library(Matrix)
library(bravo)
library(doParallel) # For parallelizing on Linux/Unix machines

# We want to run the GWAS for each environment/phenotype parallely:
# How many cores do we want to use?
# my system has 10 cores(20 threads), but we have only four environments.
# So need only 4 cores.
# Note: detectCores() give the number of available cores.
# Note: on nova then you have to request 4 CPU cores as well.
registerDoMC(cores = 1)

# Load the BLUPs
yy = read.csv("blues.csv _phenotype.csv")
str(yy)
dim(yy)
dim(Z)
# Load the numerical SNP matrix (the output of "SNP_to_X.r")
# Note it's an RDS (single data object) file hence read using readRDS.
Z = readRDS("zmat_sparse.rds")

#check if the rows of yy and rows of Z are aligned:
all.equal(rownames(Z),yy$Taxa)

# No... seems like the first three lines of BLUPs were not genotyped.
Y = yy

# check again
all.equal(rownames(Z),Y$Taxa)
# We will use SVEN, but there's a tuning parameter called "lam". The default value of 
# lam is n/p^2 where n is the number of lines and p = number of markers.
# But in certain scenarios (e.g. when markers are highly correlated) larger values of 
# lam could be useful (see Li et al., 2020)
# Run ?sven for more info


# First try with lam=0.1... it was useful for Kevin's data on mungbeans
fits = foreach(ii = 2:ncol(Y),.packages = c("bravo","Matrix")) %dopar%
  {
    X = Z + 0
    y = Y[,ii]
    phen.name = names(Y)[ii]
    set.seed(2441139)
    fit = sven(X = X,y = y,Ntemp = 30,Miter = 200,verbose = F,lam = 20)
    fit$Phenotype=phen.name
    return(fit)
  }
# save the result
saveRDS(object = fits,file="GWAS-results-1.rds")

rm(fits)


# Next move to SVEN post processing..
library(bravo)
library(dplyr)
library(tibble)
################
# We need the marker matrix to get the SNP names
Z <- readRDS("zmat_sparse.rds")
snp.names <- colnames(Z)

gwas <- readRDS("GWAS-results-1.rds")

all_results <- lapply(seq_along(gwas), function(i) {
  model <- gwas[[i]]$model.wam
  trait <- gwas[[i]]$stats$y
  Xmat  <- as.matrix(gwas[[i]]$stats$Xm[, model, drop = FALSE])
  
  # variance explained per SNP (R² * 100)
  pve <- as.numeric(cor(trait, Xmat)^2)
  pve <- round(pve * 100, 1)
  
  # MIP etc. from SVEN
  Summary <- mip.sven(gwas[[i]], threshold = -1)[model, , drop = FALSE]
  
  # move rownames (SNP IDs) into a column
  Summary <- as.data.frame(Summary)
  Summary <- rownames_to_column(Summary, var = "SNP")
  
  # add phenotype and PVE
  Summary$Phenotype <- gwas[[i]]$Phenotype
  Summary$PVE       <- pve
  
  # optional: order columns nicely
  Summary <- Summary |>
    select(Phenotype, SNP, PVE, everything())
  
  Summary
})

# bind all phenotypes into one data frame
results_all <- bind_rows(all_results)

# write to CSV
write.csv(
  results_all,
  file      = "SVEN_MIP_PVE_results_labda20_Dutta.csv",
  row.names = FALSE
)

cat("Saved SVEN summary to SVEN_MIP_PVE_results.csv\n")

#######new code_Gemini####
library(bravo)
library(dplyr)
library(tibble)
library(Matrix) 

# 1. Load Data
Z <- readRDS("zmat_sparse.rds")
gwas <- readRDS("GWAS-results-1.rds")

# 2. Run Extraction Loop
all_results <- lapply(seq_along(gwas), function(i) {
  
  # Get indices of selected SNPs
  model <- gwas[[i]]$model.wam
  
  # Safety Check: If no SNPs selected, return NULL
  if (length(model) == 0) {
    return(NULL)
  }
  
  # Get Phenotype Data
  trait <- gwas[[i]]$stats$y
  
  # A. Get Genotype Data for PVE
  # Use global Z matrix to be safe
  Xmat <- as.matrix(Z[, model, drop = FALSE])
  
  # B. Calculate PVE (Correlation^2 * 100)
  pve <- as.numeric(cor(trait, Xmat)^2)
  pve <- round(pve * 100, 2)
  
  # C. Get MIPs (Model Inclusion Probabilities)
  # threshold = 0 ensures we get stats for all selected variables
  Summary <- mip.sven(gwas[[i]], threshold = 0)
  
  # --- THE FIX IS HERE ---
  # We must subset using CHARACTER names, not integer positions
  # This matches the row name "494" instead of looking for the 494th row
  summary_subset <- Summary[as.character(model), , drop = FALSE]
  
  # D. Create Result Dataframe
  df_res <- as.data.frame(summary_subset)
  
  # Make the "Index" a column (e.g., "494")
  df_res <- rownames_to_column(df_res, var = "SNP_Index")
  
  # E. Get Actual SNP Names (e.g., "Chr1_305...") from Z matrix
  # If Z has column names, this maps the index to the name
  real_snp_names <- colnames(Z)[model]
  if(is.null(real_snp_names)) real_snp_names <- paste0("SNP_", model)
  
  df_res$SNP_Name <- real_snp_names
  df_res$Phenotype <- gwas[[i]]$Phenotype
  df_res$PVE       <- pve
  
  # Final Cleanup
  df_res <- df_res %>%
    select(Phenotype, SNP_Index, SNP_Name, MIP, PVE)
  
  return(df_res)
})

# 3. Bind and Save
results_all <- bind_rows(all_results)

if (nrow(results_all) > 0) {
  write.csv(
    results_all,
    file      = "SVEN_MIP_PVE_results_lam20_Fixed.csv",
    row.names = FALSE
  )
  cat("Success! Saved", nrow(results_all), "QTLs to 'SVEN_MIP_PVE_results_lam30_Fixed.csv'.\n")
} else {
  cat("Warning: No SNPs were selected for any trait.\n")
}
