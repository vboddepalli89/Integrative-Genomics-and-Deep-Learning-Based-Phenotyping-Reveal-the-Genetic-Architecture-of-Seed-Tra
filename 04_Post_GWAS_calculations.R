# Load necessary libraries
library(data.table)
library(dplyr)
library(tidyr)

# ==============================================================================
# 1. File Paths (Updated to your specific location)
# ==============================================================================
gwas_file   <- "C:/Users/naresh89/Desktop/My projects/gwas/Seed image analysis/manuscript/weilv9_results/GWAS_results_four_models.csv"
pheno_file  <- "C:/Users/naresh89/Desktop/My projects/gwas/Seed image analysis/manuscript/weilv9_results/blues.csv _phenotype.csv"
geno_file   <- "C:/Users/naresh89/Desktop/My projects/gwas/Seed image analysis/manuscript/weilv9_results/blues.csv _genotype.txt"

# ==============================================================================
# 2. Load Data
# ==============================================================================
cat("Loading data...\n")
gwas_results <- fread(gwas_file)
pheno_df     <- fread(pheno_file)
geno_df      <- fread(geno_file)

# ==============================================================================
# 3. Pre-processing
# ==============================================================================

# Helper function to map IUPAC Genotypes to Numeric (0, 1, 2)
# 0 = Major Homozygote, 1 = Het, 2 = Minor Homozygote
convert_geno <- function(geno_vec) {
  # Remove NAs
  valid_g <- geno_vec[geno_vec != "N" & !is.na(geno_vec)]
  
  if(length(valid_g) == 0) return(list(nums=rep(NA, length(geno_vec)), minor_allele=NA))
  
  # Count alleles
  # Simplified: counts occurrences of homozygotes to find Major/Minor
  counts <- table(valid_g)
  homozygotes <- names(counts)[names(counts) %in% c("A", "C", "G", "T")]
  
  if(length(homozygotes) < 1) return(list(nums=rep(NA, length(geno_vec)), minor_allele=NA))
  
  # Determine Major and Minor Alleles
  sorted_hom <- names(sort(counts[homozygotes], decreasing = TRUE))
  major <- sorted_hom[1]
  minor <- if(length(sorted_hom) > 1) sorted_hom[2] else NA # Handle monomorphic
  
  # Mapping rules
  # Heterozygotes map: R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C
  
  nums <- sapply(geno_vec, function(x) {
    if(is.na(x) || x == "N") return(NA)
    if(x == major) return(0)
    if(!is.na(minor) && x == minor) return(2)
    
    # Check Heterozygotes
    if(x %in% c("R","Y","S","W","K","M")) return(1)
    
    return(NA) 
  })
  
  return(list(nums=nums, minor_allele=minor))
}

# ==============================================================================
# 4. Calculation Loop
# ==============================================================================
cat("Calculating PVE and Effects...\n")

# Initialize columns
gwas_results$PVE <- NA
gwas_results$Effect <- NA
gwas_results$SE <- NA
gwas_results$Effect_Allele <- NA

# Get list of sample names from Genotype file that overlap with Phenotype file
common_taxa <- intersect(names(geno_df)[-c(1:11)], pheno_df$Taxa)

# Loop through each row in GWAS results
for(i in 1:nrow(gwas_results)) {
  
  snp_id <- gwas_results$SNP[i]
  trait  <- gwas_results$Trait[i]
  
  # 1. Get Genotype Data for this SNP
  snp_row <- geno_df[geno_df$`rs#` == snp_id, ..common_taxa]
  
  if(nrow(snp_row) == 0) next # SNP not found in genotype file
  
  geno_vec <- as.character(as.vector(t(snp_row)))
  
  # 2. Get Phenotype Data for this Trait
  pheno_vec <- pheno_df[match(common_taxa, pheno_df$Taxa)][[trait]]
  
  # 3. Convert Genotype to Numeric
  geno_data <- convert_geno(geno_vec)
  geno_num  <- geno_data$nums
  eff_allele <- geno_data$minor_allele
  
  # 4. Run Linear Model (Phenotype ~ Genotype)
  if(!all(is.na(geno_num)) && !all(is.na(pheno_vec))) {
    model <- lm(pheno_vec ~ geno_num)
    mod_sum <- summary(model)
    
    # Store Results
    gwas_results$PVE[i] <- mod_sum$r.squared
    gwas_results$Effect[i] <- coef(model)[2] # Slope (Effect of Minor Allele)
    gwas_results$SE[i] <- mod_sum$coefficients[2, "Std. Error"]
    gwas_results$Effect_Allele[i] <- eff_allele
  }
}

# ==============================================================================
# 5. Save Results
# ==============================================================================
output_file <- "C:/Users/naresh89/Desktop/My projects/gwas/Seed image analysis/manuscript/weilv9_results/GWAS_results_with_Stats.csv"
write.csv(gwas_results, output_file, row.names = FALSE)

cat("Done! Results saved to:", output_file, "\n")
