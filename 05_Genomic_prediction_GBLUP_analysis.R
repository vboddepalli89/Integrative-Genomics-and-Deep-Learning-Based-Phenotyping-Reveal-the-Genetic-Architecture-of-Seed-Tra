################################################################################
## Integrative Genomic Prediction Pipeline (gBLUP and gBLUP + SNPs)
## Evaluates genomic prediction accuracies and selection coincidence index (CI)
################################################################################

# Clear workspace
rm(list = ls())
start.time <- Sys.time()

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(rrBLUP)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(parallel)

cat("Number of available cores:", detectCores(), "\n")

# ------------------------------------------------------------------------------
# 2. Data Import and Preparation
# ------------------------------------------------------------------------------
# Load Genotype Data
geno_raw <- read.table("seed_numerical.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
cat("Original genotype dimensions:", dim(geno_raw), "\n")

# Convert to numeric matrix
X <- apply(as.matrix(geno_raw), 2, as.numeric)
rownames(X) <- rownames(geno_raw)

# Load Phenotype Data (BLUEs)
blues <- read.csv("blues.csv", stringsAsFactors = FALSE)
cat("Phenotype dimensions:", dim(blues), "\n")
cat("Phenotype traits:", colnames(blues)[-1], "\n")

# Load Trait-Specific Significant SNPs (GWAS output)
snps_all <- read.csv("sig_snps_3.csv", stringsAsFactors = FALSE)

# Match Taxa between Genotype and Phenotype
common_genos <- intersect(rownames(X), blues$Taxa)
cat("Number of common genotypes:", length(common_genos), "\n")

X <- X[common_genos, , drop = FALSE]
y <- blues[match(common_genos, blues$Taxa), ]

# ------------------------------------------------------------------------------
# 3. Compute Genomic Relationship Matrix (A Matrix)
# ------------------------------------------------------------------------------
cat("Computing genomic relationship matrix with EM imputation...\n")
A1 <- A.mat(X, impute.method = "EM", n.core = 1, max.missing = 0.5)
rownames(A1) <- common_genos
colnames(A1) <- common_genos

# ------------------------------------------------------------------------------
# 4. Clean and Align SNP List
# ------------------------------------------------------------------------------
# Remove duplicate (Trait, SNP) rows and keep only SNPs present in genotype matrix
snps_all <- snps_all[!duplicated(snps_all[, c("Trait", "SNP")]), ]
snps_all$in_geno <- snps_all$SNP %in% colnames(X)
snps_use <- snps_all[snps_all$in_geno, c("Trait", "SNP")]

# ------------------------------------------------------------------------------
# 5. Cross-Validation Setup
# ------------------------------------------------------------------------------
set.seed(123)
cycles <- 100
folds  <- 10

trait_names <- colnames(y)[-1]
n_traits <- length(trait_names)

# Storage matrices for accuracies
acc_gblup <- matrix(nrow = cycles * folds, ncol = n_traits)
acc_snp   <- matrix(nrow = cycles * folds, ncol = n_traits)
colnames(acc_gblup) <- trait_names
colnames(acc_snp)   <- trait_names

# Storage for GEBVs (first cycle)
gebv_gblup_list <- list()
gebv_snp_list   <- list()

# ------------------------------------------------------------------------------
# 6. Main Cross-Validation Loop (Evaluates Both Models Simultaneously)
# ------------------------------------------------------------------------------
cat("\nStarting", cycles, "cycles of", folds, "-fold cross-validation...\n")

for (cycle in 1:cycles) {
  if (cycle %% 10 == 0) cat(" Processing cycle", cycle, "of", cycles, "\n")
  
  set.seed(cycle)
  fold_assignments <- sample(rep(1:folds, length.out = length(common_genos)))
  
  cycle_gebv_gblup <- list()
  cycle_gebv_snp   <- list()
  
  for (trait_idx in 2:ncol(y)) {
    trait_name <- colnames(y)[trait_idx]
    y_vec <- y[[trait_name]]
    
    # Prepare Fixed Effects for gBLUP+SNPs
    trait_snps <- snps_use$SNP[snps_use$Trait == trait_name]
    if (length(trait_snps) > 0) {
      X_fixed_base <- cbind(Intercept = 1, X[common_genos, trait_snps, drop = FALSE])
    } else {
      X_fixed_base <- cbind(Intercept = 1)
    }
    
    for (fold in 1:folds) {
      test_genos <- common_genos[fold_assignments == fold]
      y_test_actual <- y_vec[match(test_genos, common_genos)]
      row_idx <- (cycle - 1) * folds + fold
      
      # --- MODEL 1: Standard gBLUP ---
      y_train_df <- y
      y_train_df[y_train_df$Taxa %in% test_genos, trait_idx] <- NA
      
      fit_gblup <- kin.blup(data = y_train_df, geno = "Taxa", pheno = trait_name, K = A1)
      y_test_pred_gblup <- fit_gblup$g[test_genos]
      acc_gblup[row_idx, trait_idx - 1] <- cor(y_test_pred_gblup, y_test_actual, use = "complete.obs")
      
      if (cycle == 1 && fold == 1) cycle_gebv_gblup[[trait_name]] <- fit_gblup$g
      
      # --- MODEL 2: gBLUP + SNPs ---
      y_train_vec <- y_vec
      y_train_vec[match(test_genos, common_genos)] <- NA
      
      fit_snp <- mixed.solve(y = y_train_vec, K = A1, X = X_fixed_base)
      u_hat <- fit_snp$u
      names(u_hat) <- common_genos
      y_fixed_hat <- as.vector(X_fixed_base %*% fit_snp$beta)
      
      y_test_pred_snp <- y_fixed_hat[match(test_genos, common_genos)] + u_hat[test_genos]
      acc_snp[row_idx, trait_idx - 1] <- cor(y_test_pred_snp, y_test_actual, use = "complete.obs")
      
      if (cycle == 1 && fold == 1) {
        trait_gebv_snp <- y_fixed_hat + u_hat
        names(trait_gebv_snp) <- common_genos
        cycle_gebv_snp[[trait_name]] <- trait_gebv_snp
      }
    }
  }
  
  if (cycle == 1) {
    gebv_gblup_list <- cycle_gebv_gblup
    gebv_snp_list   <- cycle_gebv_snp
  }
}

# ------------------------------------------------------------------------------
# 7. Summarize Prediction Accuracies
# ------------------------------------------------------------------------------
sum_gblup <- data.frame(
  Trait = colnames(acc_gblup),
  Mean_Acc_gBLUP = apply(acc_gblup, 2, mean, na.rm = TRUE),
  SD_Acc_gBLUP   = apply(acc_gblup, 2, sd, na.rm = TRUE)
)

sum_snp <- data.frame(
  Trait = colnames(acc_snp),
  Mean_Acc_SNP = apply(acc_snp, 2, mean, na.rm = TRUE),
  SD_Acc_SNP   = apply(acc_snp, 2, sd, na.rm = TRUE)
)

comp_tbl <- merge(sum_gblup, sum_snp, by = "Trait") %>%
  mutate(Diff_in_Accuracy = Mean_Acc_SNP - Mean_Acc_gBLUP) %>%
  arrange(desc(Diff_in_Accuracy))

cat("\n--- Accuracy Comparison Summary ---\n")
print(comp_tbl)

# ------------------------------------------------------------------------------
# 8. Selection Coincidence Index (Top 10 Genotypes)
# ------------------------------------------------------------------------------
get_topk <- function(x, taxa, k = 10, decreasing = TRUE) {
  ord <- order(x, decreasing = decreasing, na.last = NA)
  taxa[ord][seq_len(min(k, sum(!is.na(x))))]
}

ci_list <- list()
taxa_vec <- y$Taxa

for (tr in trait_names) {
  blup_vals  <- y[[tr]]
  top_blup   <- get_topk(blup_vals, taxa = taxa_vec, k = 10)
  
  # Standard gBLUP CI
  top_gblup  <- get_topk(gebv_gblup_list[[tr]], taxa = taxa_vec, k = 10)
  ci_gblup   <- length(intersect(top_blup, top_gblup)) / length(top_blup) * 100
  
  # gBLUP+SNP CI
  top_snp    <- get_topk(gebv_snp_list[[tr]], taxa = taxa_vec, k = 10)
  ci_snp_val <- length(intersect(top_blup, top_snp)) / length(top_blup) * 100
  
  ci_list[[tr]] <- data.frame(
    Trait = tr,
    GBLUP_CI = ci_gblup,
    GBLUP_SNP_CI = ci_snp_val
  )
}

ci_results <- bind_rows(ci_list)

# ------------------------------------------------------------------------------
# 9. Plotting Results (ggplot2 & patchwork)
# ------------------------------------------------------------------------------
trait_labels <- c(
  "Seed_area_mm2"     = "Seed area",
  "Seed_length_mm"    = "Seed length",
  "Seed_width_mm"     = "Seed width",
  "Aspect_Ratio"      = "Seed aspect ratio",
  "Seed_perimeter_mm" = "Seed perimeter",
  "Seed_circularity"  = "Seed circularity",
  "SWT100"            = "100 Seed weight"
)

# Prepare Accuracy Data for Plotting
acc_long_gblup <- data.frame(Trait = rep(colnames(acc_gblup), each = nrow(acc_gblup)), Accuracy = as.vector(acc_gblup), Model = "GBLUP")
acc_long_snp   <- data.frame(Trait = rep(colnames(acc_snp), each = nrow(acc_snp)), Accuracy = as.vector(acc_snp), Model = "GBLUP+SNP")
acc_combined   <- rbind(acc_long_gblup, acc_long_snp)

p_acc <- ggplot(acc_combined, aes(x = Trait, y = Accuracy, fill = Model)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  coord_flip() +
  labs(title = "A. Genomic Prediction Accuracy: GBLUP vs GBLUP+SNP", x = "Trait", y = "Prediction Accuracy") +
  scale_fill_manual(values = c("GBLUP" = "grey70", "GBLUP+SNP" = "steelblue")) +
  scale_x_discrete(labels = trait_labels) +
  theme_bw()

# Prepare CI Data for Plotting
ci_long <- ci_results %>%
  pivot_longer(cols = c("GBLUP_CI", "GBLUP_SNP_CI"), names_to = "Model", values_to = "CI") %>%
  mutate(Model = ifelse(Model == "GBLUP_CI", "GBLUP", "GBLUP+SNP"))

p_ci <- ggplot(ci_long, aes(x = reorder(Trait, CI), y = CI, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7)) +
  labs(title = "B. Selection Coincidence Index (Top 10)", x = "Trait", y = "Coincidence Index (%)") +
  scale_fill_manual(values = c("GBLUP" = "grey70", "GBLUP+SNP" = "steelblue")) +
  scale_x_discrete(labels = trait_labels) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine and Save Plot
combined_plot <- p_acc / p_ci
print(combined_plot)
ggsave("Combined_GP_Comparison_Plot.jpg", plot = combined_plot, width = 10, height = 12, dpi = 300)

# ------------------------------------------------------------------------------
# 10. Save CSV Exports
# ------------------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M")
write.csv(comp_tbl, paste0("GP_Accuracy_Summary_Comparison_", timestamp, ".csv"), row.names = FALSE)
write.csv(ci_results, paste0("GP_Coincidence_Index_Comparison_", timestamp, ".csv"), row.names = FALSE)

cat("\nPipeline complete! Total time:", round(difftime(Sys.time(), start.time, units = "mins"), 2), "minutes\n")