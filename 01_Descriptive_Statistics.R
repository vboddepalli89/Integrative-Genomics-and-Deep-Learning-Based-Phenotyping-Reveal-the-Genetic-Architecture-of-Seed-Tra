# ==============================================================================
# Seed Image Analysis - GWAS BLUEs Generation
# ==============================================================================

# 1. Setup and Libraries -------------------------------------------------------
setwd("C:/Users/naresh89/Desktop/My projects/gwas/Seed image analysis/descriptive stats")

# Check if packages are installed, if not, install them
required_packages <- c("lme4", "lmerTest", "emmeans", "readr", "dplyr", "car", "MASS", "openxlsx")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(lme4)
library(lmerTest)
library(emmeans)
library(readr)
library(dplyr)
library(car)
library(MASS)
library(openxlsx)

# 2. Load Data -----------------------------------------------------------------
data <- read.xlsx("seed_traits_2022_24_final.xlsx")

# Ensure factors are set correctly
data$line <- as.factor(data$line)
data$env <- as.factor(data$env)
data$block <- as.factor(data$block)

# List of traits to analyze
traits <- c("Seed_area_mm2", "Seed_length_mm", "Seed_width_mm", 
            "Aspect_Ratio", "Seed_perimeter_mm", "Seed_circularity")

# Initialize a list to store results
blues_list <- list()

# 3. Analysis Loop (OPTIMIZED) -------------------------------------------------

# Set emmeans options to use asymptotic method (FAST)
# This prevents the code from getting stuck calculating degrees of freedom
emm_options(pbkrtest.limit = 0, lmerTest.limit = 0)

for (trait in traits) {
  cat(paste0("\nProcessing Trait: ", trait, "...\n"))
  
  # Create a temporary dataframe for specific trait analysis (remove NAs)
  df_trait <- data %>% 
    dplyr::select(line, env, block, all_of(trait)) %>% 
    na.omit()
  
  colnames(df_trait)[4] <- "value"
  
  # --- A. Outlier Removal ---
  outlier_model <- lm(value ~ line + env + block, data = df_trait)
  outlier_test <- outlierTest(outlier_model, cutoff = 0.05, n.max = 10)
  
  if (!is.null(outlier_test$rstudent)) {
    outliers <- names(outlier_test$rstudent)
    if (length(outliers) > 0) {
      cat(paste0("  - Found ", length(outliers), " outliers. Removing...\n"))
      df_trait <- df_trait[-as.numeric(outliers), ]
    }
  } else {
    cat("  - No significant outliers detected.\n")
  }
  
  # --- B. Normality Check & Transformation ---
  
  # >>>>> START OF NEW CODE <<<<<
  # Define traits that must remain raw (Physical dimensions)
  # You can add "Aspect_Ratio" here if you prefer raw ratios over transformed ones
  physical_traits <- c("Seed_length_mm", "Seed_width_mm", "Seed_perimeter_mm", "Seed_area_mm2")
  
  if (trait %in% physical_traits) {
    cat(paste0("  - Physical dimension detected ('", trait, "'). FORCE KEEPING RAW DATA.\n"))
    transformation_applied <- "None"
    
  } else {
    # >>>>> END OF NEW CODE (Start of 'ELSE' block) <<<<<
    
    # --- EXISTING NORMALITY LOGIC GOES INSIDE THIS ELSE BLOCK ---
    
    # Use a simple model for normality check to speed up
    lm_norm <- lm(value ~ line + env, data = df_trait)
    residuals <- resid(lm_norm)
    
    # Sample residuals if too many
    if(length(residuals) > 5000) residuals <- sample(residuals, 5000)
    shapiro_res <- shapiro.test(residuals)
    
    transformation_applied <- "None"
    
    if (shapiro_res$p.value < 0.05) {
      cat(paste0("  - Normality test significant (p < 0.05). Checking Box-Cox...\n"))
      
      bc <- boxcox(lm(value ~ 1, data = df_trait), plotit = FALSE)
      lambda <- bc$x[which.max(bc$y)]
      
      cat(paste0("  - Optimal Lambda: ", round(lambda, 2), "\n"))
      
      if (lambda > -0.2 & lambda < 0.2) {
        df_trait$value <- log(df_trait$value)
        transformation_applied <- "Log"
      } else if (lambda > 0.3 & lambda < 0.7) {
        df_trait$value <- sqrt(df_trait$value)
        transformation_applied <- "Sqrt"
      } else if (lambda < -0.8) {
        df_trait$value <- 1 / df_trait$value
        transformation_applied <- "Inverse"
      } else {
        cat("  - Lambda close to 1 or complex. Keeping raw data.\n")
      }
      cat(paste0("  - Applied Transformation: ", transformation_applied, "\n"))
    }
    
  } 
  
  # --- C. Generate BLUEs (Mixed Model) ---
  cat("  - Fitting Mixed Model (this may take a moment)...\n")
  
  # Using lme4::lmer explicitly. 
  # Added control optimizer to speed up convergence
  final_model <- lme4::lmer(value ~ line + (1|env) + (1|env:block) + (1|line:env), 
                            data = df_trait,
                            control = lmerControl(optimizer = "bobyqa"))
  
  cat("  - Calculating BLUEs...\n")
  # With limits set to 0 above, this should now be instant
  blues <- emmeans(final_model, specs = "line")
  
  blues_df <- as.data.frame(blues)
  
  blues_clean <- blues_df %>% 
    dplyr::select(line, emmean) %>% 
    rename(!!trait := emmean)
  
  # Store in list
  blues_list[[trait]] <- blues_clean
}

# ==============================================================================
# Generate Data Processing Summary Table
# ==============================================================================

# Initialize a data frame to store the stats
processing_summary <- data.frame(
  Trait = character(),
  Outliers_Removed = integer(),
  Shapiro_P_Value = numeric(),
  BoxCox_Lambda = numeric(),
  Transformation_Used = character(),
  stringsAsFactors = FALSE
)

cat("\nGenerating Data Processing Summary...\n")

for (trait in traits) {
  
  # 1. Prepare Data (Same as before)
  df_trait <- data %>% 
    dplyr::select(line, env, block, all_of(trait)) %>% 
    na.omit()
  colnames(df_trait)[4] <- "value"
  
  # 2. Count Outliers
  outlier_model <- lm(value ~ line + env + block, data = df_trait)
  outlier_test <- outlierTest(outlier_model, cutoff = 0.05, n.max = 10)
  
  n_outliers <- 0
  if (!is.null(outlier_test$rstudent)) {
    outliers <- names(outlier_test$rstudent)
    n_outliers <- length(outliers)
    if (n_outliers > 0) {
      df_trait <- df_trait[-as.numeric(outliers), ] # Remove for normality check
    }
  }
  
  # 3. Check Normality & Lambda
  lm_norm <- lm(value ~ line + env, data = df_trait)
  residuals <- resid(lm_norm)
  if(length(residuals) > 5000) residuals <- sample(residuals, 5000)
  shapiro_res <- shapiro.test(residuals)
  
  # Calculate Lambda (Box-Cox) regardless of p-value for reporting
  bc <- boxcox(lm(value ~ 1, data = df_trait), plotit = FALSE)
  lambda <- bc$x[which.max(bc$y)]
  
  # Determine what transformation WAS applied in your main loop
  trans_applied <- "None"
  if (shapiro_res$p.value < 0.05) {
    if (lambda > -0.2 & lambda < 0.2) {
      trans_applied <- "Log"
    } else if (lambda > 0.3 & lambda < 0.7) {
      trans_applied <- "Sqrt"
    } else if (lambda < -0.8) {
      trans_applied <- "Inverse"
    }
  }
  
  # 4. Add to Summary Table
  processing_summary[nrow(processing_summary) + 1, ] <- list(
    Trait = trait,
    Outliers_Removed = n_outliers,
    Shapiro_P_Value = signif(shapiro_res$p.value, 3), # Round for readability
    BoxCox_Lambda = round(lambda, 2),
    Transformation_Used = trans_applied
  )
}

# Print the table to console
print(processing_summary)

# Save to CSV
write.csv(processing_summary, "Data_Processing_Summary.csv", row.names = FALSE)

cat("\nSummary saved to 'Data_Processing_Summary.csv'\n")
# 4. Merge and Save ------------------------------------------------------------
cat("\nMerging all traits...\n")

# Merge all dataframes in the list by 'line'
final_blues <- Reduce(function(x, y) merge(x, y, by = "line", all = TRUE), blues_list)

# View first few rows
print(head(final_blues))

# Save to CSV
output_filename <- "Seed_Size_GWAS_BLUEs_Feb.csv"
write.csv(final_blues, output_filename, row.names = FALSE)

cat(paste0("\nAnalysis Complete. BLUEs saved to: ", output_filename, "\n"))
# ==============================================================================
# 5. Descriptive Statistics Generation
# ==============================================================================

# Ensure tidyr is loaded for data reshaping
if(!"tidyr" %in% installed.packages()[,"Package"]) install.packages("tidyr")
library(tidyr)

cat("\nCalculating descriptive statistics...\n")

# 1. Calculate Statistics
# We pivot the data to 'long' format to calculate stats for all traits at once
desc_stats <- final_blues %>%
  pivot_longer(cols = -line, names_to = "Trait", values_to = "BLUE_Value") %>%
  group_by(Trait) %>%
  summarise(
    N_Lines = n(),
    Mean = mean(BLUE_Value, na.rm = TRUE),
    Min = min(BLUE_Value, na.rm = TRUE),
    Max = max(BLUE_Value, na.rm = TRUE),
    SD = sd(BLUE_Value, na.rm = TRUE),
    CV_Percent = (sd(BLUE_Value, na.rm = TRUE) / mean(BLUE_Value, na.rm = TRUE)) * 100
  ) %>%
  mutate_if(is.numeric, round, 2) # Round all numeric columns to 2 decimal places

print(desc_stats)

# 2. Save to Excel (Professional Format)
# Since you loaded 'openxlsx', let's use it to save both the BLUEs and Stats 
# in a single file on different sheets (very convenient for publication).

wb <- createWorkbook()

# Add Sheet 1: The BLUEs (Raw Data)
addWorksheet(wb, "Genotype BLUEs")
writeData(wb, "Genotype BLUEs", final_blues)

# Add Sheet 2: The Descriptive Statistics (Summary)
addWorksheet(wb, "Descriptive Stats")
writeData(wb, "Descriptive Stats", desc_stats)

# Save the workbook
saveWorkbook(wb, "Seed_Size_GWAS_Results_Feb_Complete.xlsx", overwrite = TRUE)

cat(paste0("\nSuccess! Stats and BLUEs saved to: Seed_Size_GWAS_Results_Complete.xlsx\n"))
# ==============================================================================
# 6. Heritability & Variance Component Analysis
# ==============================================================================

# Initialize a list to store heritability results
h2_results_list <- list()

cat("\nStarting Heritability Analysis (with Outlier Removal)...\n")

for (trait in traits) {
  
  # 1. Prepare Data
  df_trait <- data %>% 
    dplyr::select(line, env, block, all_of(trait)) %>% 
    na.omit()
  colnames(df_trait)[4] <- "value"
  
  # --- ADDED: Outlier Removal (Same as BLUEs loop) ---
  outlier_model <- lm(value ~ line + env + block, data = df_trait)
  outlier_test <- car::outlierTest(outlier_model, cutoff = 0.05, n.max = 10)
  
  if (!is.null(outlier_test$rstudent)) {
    outliers <- names(outlier_test$rstudent)
    if (length(outliers) > 0) {
      df_trait <- df_trait[-as.numeric(outliers), ]
    }
  }
  # ---------------------------------------------------
  
  # 2. Fit Random Effects Model
  if (n_distinct(df_trait$env) > 1) {
    model_h2 <- lmer(value ~ (1|line) + (1|env) + (1|line:env) + (1|env:block), data = df_trait)
  } else {
    model_h2 <- lmer(value ~ (1|line) + (1|block), data = df_trait)
  }
  
  # 3. Extract Variance Components & Calculate H2
  vc <- as.data.frame(VarCorr(model_h2))
  sigma_g   <- vc[vc$grp == "line", "vcov"]
  sigma_e   <- vc[vc$grp == "Residual", "vcov"]
  
  if (n_distinct(df_trait$env) > 1) {
    sigma_ge  <- vc[vc$grp == "line:env", "vcov"]
    if(length(sigma_ge)==0) sigma_ge <- 0
  } else {
    sigma_ge <- 0
  }
  
  n_env <- n_distinct(df_trait$env)
  reps_per_line <- table(df_trait$line)
  n_rep <- length(reps_per_line) / sum(1/reps_per_line)
  
  denominator <- sigma_g + (sigma_ge / n_env) + (sigma_e / (n_env * n_rep))
  H2 <- sigma_g / denominator
  
  # Store result
  h2_results_list[[trait]] <- data.frame(
    Trait = trait,
    Vg = round(sigma_g, 4),
    Vge = round(sigma_ge, 4),
    Ve = round(sigma_e, 4),
    H2_BroadSense = round(H2, 2)
  )
}

# Print Updated Table
h2_table_corrected <- do.call(rbind, h2_results_list)
print(h2_table_corrected)

# 7. Add to your Excel file (Append to existing workbook)
# We load the workbook we created in Step 5 and add this new sheet
wb <- loadWorkbook("Seed_Size_GWAS_Results_Feb_Complete.xlsx")
addWorksheet(wb, "Heritability & Variances")
# Use the correct variable name: h2_table_corrected
writeData(wb, "Heritability & Variances", h2_table_corrected)

# Save the workbook
saveWorkbook(wb, "Seed_Size_GWAS_Results_Feb_Complete.xlsx", overwrite = TRUE)

cat("\n✅ Success! All results (BLUEs, Stats, Heritability) are saved.\n")
