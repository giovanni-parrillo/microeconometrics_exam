library(fixest)
library(lubridate)
library(dplyr)

# Load the fixed dataset
# Make sure this is the file you just created with the post-years included
matched_data <- read.csv("matched_data_fixed.csv")
panel_data <- read.csv("rta_panel.csv")
out_file <- "tabelle_latex/regression_tables.tex"

# Initialize the LaTeX file with preamble
cat("\\documentclass{article}
\\usepackage{booktabs}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{geometry}
\\usepackage{graphicx}
\\geometry{a4paper, margin=1in}

\\begin{document}

", file = out_file)


# Create Variables
# Create Post variable (1 if year >= entry_year)
matched_data$post_rta <- ifelse(matched_data$year >= matched_data$Entry.into.Force, 1, 0)
panel_data$post_rta <- ifelse(panel_data$year >= panel_data$Entry.into.Force, 1, 0) # nolint
# Create Log Loss (Outcome)
matched_data$ln_loss <- log(matched_data$loss+1)
panel_data$ln_loss <- log(panel_data$loss+1) #because in panel_data there are losses=0

# Diagnostic Check
print("--- Distribution check ---")
print(table(Enviro = matched_data$enviro_rta, Post = matched_data$post_rta))

print(table(Enviro = panel_data$enviro_rta, Post = panel_data$post_rta))

#Robust covariance matrix with sparse cross-cluster correlation

w <- function(g, h) {
  # Filter matched_data for the given id g
  g_data <- matched_data[matched_data$id == g, ][1, 2:202]
  h_data <- matched_data[matched_data$id == h, ][1, 2:202]

  # Sum the number of 1s from column 2 to column 202 for id g
  sum_g <- sum(g_data, na.rm = TRUE)
  # Count columns where both g and h have 1 (between columns 2 and 202)
  both_ones <- sum((g_data == 1) & (h_data == 1), na.rm = TRUE)
  # Divide the count of common 1s by the sum of g's 1s
  result <- both_ones / sum_g
  
  return(result)
}

w_gen <- function(dataset, g, h) {
  # Filter matched_data for the given id g
  g_data <- dataset[dataset$id == g, ][1, 2:202]
  h_data <- dataset[dataset$id == h, ][1, 2:202]
  
  # Sum the number of 1s from column 2 to column 202 for id g
  sum_g <- sum(g_data, na.rm = TRUE)
  # Count columns where both g and h have 1 (between columns 2 and 202)
  both_ones <- sum((g_data == 1) & (h_data == 1), na.rm = TRUE)
  # Divide the count of common 1s by the sum of g's 1s
  result <- both_ones / sum_g
  
  return(result)
}

unique_ids <- unique(matched_data$id)


print("--- Generating Weight Matrix for Matched Data ---")
W_matrix <- matrix(0, length(unique_ids), length(unique_ids))

for (i in 1:length(unique_ids)){
  for (k in 1:length(unique_ids)){
    W_matrix[i,k] <- w_gen(matched_data, unique_ids[i],unique_ids[k])
  }
}


unique_ids_panel <- unique(panel_data$id)

print("--- Generating Weight Matrix for Full Panel Data ---")
W_matrix_not_matched <- matrix(0, length(unique_ids_panel), length(unique_ids_panel))

for (i in 1:length(unique_ids_panel)){
  for (k in 1:length(unique_ids_panel)){
    W_matrix_not_matched[i,k] <- w_gen(panel_data, unique_ids_panel[i],unique_ids_panel[k])
  }
}


# Function to compute robust weighted covariance matrix
# This function takes a fitted model and weight matrix to compute
# a robust covariance matrix accounting for cross-cluster correlation
compute_robust_vcov <- function(model, W_matrix, data, unique_ids) {
  # Extract residuals and design matrix from the model
  eps <- residuals(model)
  X <- model.matrix(model)
  
  # Compute observation-level scores (residuals * regressors)
  obs_scores <- sweep(X, 1, eps, '*')
  
  # Aggregate scores by cluster/group
  group_indices <- match(data$id, unique_ids)
  S_matrix <- rowsum(obs_scores, group = group_indices)
  
  # Compute the middle part of the sandwich estimator
  V_hat <- t(S_matrix) %*% W_matrix %*% S_matrix
  
  # Compute the full sandwich covariance matrix: (X'X)^-1 * V_hat * (X'X)^-1
  V <- solve(t(X) %*% X) %*% V_hat %*% solve(t(X) %*% X)
  
  return(V)
}

# Run Model to get residuals
main_ln_loss <- feols(ln_loss ~ post_rta + post_rta:enviro_rta | id + year,
               data = matched_data,
               cluster = ~id,
               weights = ~weights)

main_ln_loss_not_matched <- feols(ln_loss ~ post_rta + post_rta:enviro_rta | id + year,
                           data = panel_data,
                           cluster = ~id)

main_rate <- feols(avg_loss_rate ~ post_rta + post_rta:enviro_rta | id + year,
                      data = matched_data,
                      cluster = ~id,
                      weights = ~weights)

main_rate_not_matched <- feols(avg_loss_rate ~ post_rta + post_rta:enviro_rta | id + year,
                                  data = panel_data,
                                  cluster = ~id)

# Compute robust covariance matrix for  matched and not matched ln_loss and rate models
V_main_ln_loss <- compute_robust_vcov(main_ln_loss, W_matrix, matched_data, unique_ids)
V_main_ln_loss_not_matched <- compute_robust_vcov(main_ln_loss_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)
V_main_rate <- compute_robust_vcov(main_rate, W_matrix, matched_data, unique_ids)
V_main_rate_not_matched <- compute_robust_vcov(main_rate_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)


#View Results
# Table 2: Forest Loss and Average Loss Rate Comparison
cat("\\section*{Table 2: Forest Loss and Average Loss Rate Comparison}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

      etable(
        main_ln_loss,
        main_ln_loss_not_matched,
        main_rate,
        main_rate_not_matched,
        signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
        vcov = list(
          V_main_ln_loss, 
          V_main_ln_loss_not_matched, 
          V_main_rate, 
          V_main_rate_not_matched
        ),
        dict = c(
          ln_loss = "Log Forest Loss",
          avg_loss_rate = "Average Loss Rate",
          post_rta = "Standard RTA Effect",
          "post_rta:enviro_rta" = "Enviro Provision Effect"
        ),
        
        fixef_sizes = FALSE,
        
        tex = TRUE,
        file = out_file,
        replace = FALSE,
        title = NULL
      )

cat("}

\\vspace{2em}

", file = out_file, append = TRUE)

#Other meaningful regressions



############################################################
  # Forest loss effects by country 

  # Log Tropical Loss
  matched_data$ln_tropical_loss <- log(matched_data$tropical_loss + 1)
  panel_data$ln_tropical_loss <- log(panel_data$tropical_loss + 1)


  tropical <- feols(ln_tropical_loss ~ post_rta + post_rta:enviro_rta | id + year, 
        data = matched_data,
        cluster = ~id,
        weights = ~weights)

    V_tropical_ln_loss <- compute_robust_vcov(tropical, W_matrix, matched_data, unique_ids)


  tropical_not_matched <- feols(ln_tropical_loss ~ post_rta + post_rta:enviro_rta | id + year,
          data = panel_data,
          cluster = ~id)

          V_tropical_ln_loss_not_matched <- compute_robust_vcov(tropical_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)




  # Log Development Loss
  matched_data$ln_dev_loss <- log(matched_data$dev_loss + 1)
  panel_data$ln_dev_loss <- log(panel_data$dev_loss + 1)

  dev <- feols(ln_dev_loss ~ post_rta + post_rta:enviro_rta | id + year, 
        data = matched_data,
        cluster = ~id,
        weights = ~weights)


    V_dev_ln_loss <- compute_robust_vcov(dev, W_matrix, matched_data, unique_ids)

  dev_not_matched <- feols(ln_dev_loss ~ post_rta + post_rta:enviro_rta | id + year,
          data = panel_data,
          cluster = ~id)
            V_dev_ln_loss_not_matched <- compute_robust_vcov(dev_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)   

  # Log High Biodiversity Loss
  matched_data$ln_high_biodiversity_loss <- log(matched_data$hi_biodiv_loss + 1)
  panel_data$ln_high_biodiversity_loss <- log(panel_data$hi_biodiv_loss + 1)

  biodiversity <- feols(ln_high_biodiversity_loss ~ post_rta + post_rta:enviro_rta | id + year, 
        data = matched_data,
        cluster = ~id,
        weights = ~weights)

    V_biodiversity_ln_loss <- compute_robust_vcov(biodiversity, W_matrix, matched_data, unique_ids)

  biodiversity_not_matched <- feols(ln_high_biodiversity_loss ~ post_rta + post_rta:enviro_rta | id + year,
            data = panel_data,
            cluster = ~id)

            V_biodiversity_ln_loss_not_matched <- compute_robust_vcov(biodiversity_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)


############################################################
# Table 3A: Heterogeneity Analysis - Tropical, Developing, High Biodiversity
cat("\\section*{Table 3A: Heterogeneity Analysis - Tropical, Developing, High Biodiversity}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

# Create LaTeX summary table for tropical, developing, and high biodiversity losses
etable(
  tropical, 
  dev, 
  biodiversity,
  signif.code = c("***"=0.01, "**"=0.05, "*"=0.10), # Corrected argument name
  vcov = list(V_tropical_ln_loss, V_dev_ln_loss, V_biodiversity_ln_loss),
  dict = c(
    ln_tropical_loss = "Log Tropical Forest Loss",
    ln_dev_loss = "Log Developing Forest Loss",
    ln_high_biodiversity_loss = "Log High Biodiversity Forest Loss",
    post_rta = "Standard RTA Effect",
    "post_rta:enviro_rta" = "Enviro Provision Effect"
  ),
  
  fixef_sizes = FALSE,
  
  tex = TRUE,
  file = out_file,
  replace = FALSE,
  title = NULL
)

cat("}

\\vspace{2em}

", file = out_file, append = TRUE)


############################################################


# Log Non-Tropical Loss
matched_data$ln_non_tropical_loss <- log(matched_data$loss - matched_data$tropical_loss + 1)
panel_data$ln_non_tropical_loss <- log(panel_data$loss - panel_data$tropical_loss + 1)

non_tropical <- feols(ln_non_tropical_loss ~ post_rta + post_rta:enviro_rta | id + year, 
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_non_tropical_ln_loss <- compute_robust_vcov(non_tropical, W_matrix, matched_data, unique_ids)

non_tropical_not_matched <- feols(ln_non_tropical_loss ~ post_rta + post_rta:enviro_rta | id + year,
    data = panel_data,
    cluster = ~id)

    V_non_tropical_ln_loss_not_matched <- compute_robust_vcov(non_tropical_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)

# Log Developed Loss
matched_data$ln_developed_loss <- log(matched_data$loss - matched_data$dev_loss + 1)
panel_data$ln_developed_loss <- log(panel_data$loss - panel_data$dev_loss + 1)


developed <- feols(ln_developed_loss ~ post_rta + post_rta:enviro_rta | id + year, 
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_developed_ln_loss <- compute_robust_vcov(developed, W_matrix, matched_data, unique_ids)   

developed_not_matched <- feols(ln_developed_loss ~ post_rta + post_rta:enviro_rta | id + year,
    data = panel_data,
    cluster = ~id)

    V_developed_ln_loss_not_matched <- compute_robust_vcov(developed_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)

# Log Lower Biodiversity Loss
matched_data$ln_lower_biodiversity_loss <- log(matched_data$loss - matched_data$hi_biodiv_loss + 1)
panel_data$ln_lower_biodiversity_loss <- log(panel_data$loss - panel_data$hi_biodiv_loss + 1)

lower_biodiversity <- feols(ln_lower_biodiversity_loss ~ post_rta + post_rta:enviro_rta | id + year, 
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_lower_biodiversity_ln_loss <- compute_robust_vcov(lower_biodiversity, W_matrix, matched_data, unique_ids)

lower_biodiversity_not_matched <- feols(ln_lower_biodiversity_loss ~ post_rta + post_rta:enviro_rta | id + year,
      data = panel_data,
      cluster = ~id)

        V_lower_biodiversity_ln_loss_not_matched <- compute_robust_vcov(lower_biodiversity_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)


############################################################
# Table 3B: Heterogeneity Analysis - Non-Tropical, Developed, Lower Biodiversity
cat("\\section*{Table 3B: Heterogeneity - Non-Tropical, Developed, Lower Biodiversity}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

      # Create LaTeX summary table for non-tropical, developed, and lower biodiversity losses
    etable(
      non_tropical, 
      developed, 
      lower_biodiversity,
      signif.code = c("***"=0.01, "**"=0.05, "*"=0.10),
      vcov = list(V_non_tropical_ln_loss, V_developed_ln_loss, V_lower_biodiversity_ln_loss),
      dict = c(
        ln_non_tropical_loss = "Log Non-Tropical Forest Loss",
        ln_developed_loss = "Log Developed Forest Loss",
        ln_lower_biodiversity_loss = "Log Lower Biodiversity Forest Loss",
        post_rta = "Standard RTA Effect",
        "post_rta:enviro_rta" = "Enviro Provision Effect"
      ),
      tex = TRUE,
      file = out_file,
      replace = FALSE,
      title = NULL
    )

    cat("}

\\vspace{2em}

", file = out_file, append = TRUE)


############################################################
# Regressions on developing "dev_harvest_ha","dev_harvest_ton","dev_harvest_yield","dev_ag_exp_val", "dev_ag_exp_val"/"dev_ag_exp_ton","dev_for_prod_output", "dev_for_prod_exports"

# Regressions on developing country variables (matched data only)
matched_data$ln_dev_harvest_ha <- log(matched_data$dev_harvest_ha + 1)
matched_data$ln_dev_harvest_ton <- log(matched_data$dev_harvest_ton + 1)
matched_data$ln_dev_harvest_yield <- log(matched_data$dev_harvest_yield + 1)
matched_data$ln_dev_ag_exp_val <- log(matched_data$dev_ag_exp_val + 1)
matched_data$ln_dev_ag_exp_unit_val <- log(
  matched_data$dev_ag_exp_val / matched_data$dev_ag_exp_ton + 1)
matched_data[is.na(matched_data)] <- 0
matched_data$ln_dev_for_prod_exports <- log(
  matched_data$dev_for_prod_exports + 1)
  matched_data$ln_dev_for_prod_output <- log(
    matched_data$dev_for_prod_output + 1)


dev_harvest_ha_model <- feols(
  ln_dev_harvest_ha ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_dev_harvest_ha <- compute_robust_vcov(dev_harvest_ha_model, W_matrix, matched_data, unique_ids)

dev_harvest_ton_model <- feols(
  ln_dev_harvest_ton ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_dev_harvest_ton <- compute_robust_vcov(dev_harvest_ton_model, W_matrix, matched_data, unique_ids)

dev_harvest_yield_model <- feols(
  ln_dev_harvest_yield ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_dev_harvest_yield <- compute_robust_vcov(dev_harvest_yield_model, W_matrix, matched_data, unique_ids)

dev_ag_exp_val_model <- feols(
  ln_dev_ag_exp_val ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_dev_ag_exp_val <- compute_robust_vcov(dev_ag_exp_val_model, W_matrix, matched_data, unique_ids)

dev_ag_exp_unit_val_model <- feols(
  ln_dev_ag_exp_unit_val ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_dev_ag_exp_unit_val <- compute_robust_vcov(dev_ag_exp_unit_val_model, W_matrix, matched_data, unique_ids)

dev_for_prod_outputs_model <- feols(
  ln_dev_for_prod_output ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_dev_for_prod_output <- compute_robust_vcov(dev_for_prod_outputs_model, W_matrix, matched_data, unique_ids)

dev_for_prod_exports_model <- feols(
  ln_dev_for_prod_exports ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_dev_for_prod_exports <- compute_robust_vcov(dev_for_prod_exports_model, W_matrix, matched_data, unique_ids)

################################################
# Table 4a: Mechanism Analysis - Developing Countries
cat("\\section*{Table 4A: Mechanism - Developing Countries Agriculture \\& Forestry}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

etable(
  dev_harvest_ha_model,
  dev_harvest_ton_model,
  dev_harvest_yield_model,
  dev_ag_exp_val_model,
  dev_ag_exp_unit_val_model,
  dev_for_prod_outputs_model,
  dev_for_prod_exports_model,
  signif.code = c("***"=0.01, "**"=0.05, "*"=0.10),
  vcov = list(
    V_dev_harvest_ha,
    V_dev_harvest_ton,
    V_dev_harvest_yield,
    V_dev_ag_exp_val,
    V_dev_ag_exp_unit_val,
    V_dev_for_prod_output,
    V_dev_for_prod_exports
  ),
  dict = c(
    dev_harvest_ha = "Harvest Area (ha)",
    dev_harvest_ton = "Harvest (tons)",
    dev_harvest_yield = "Harvest Yield",
    dev_ag_exp_val = "Ag Export Value",
    dev_ag_exp_unit_val = "Ag Export Unit Value",
    dev_for_prod_output = "Forest Product Output",
    dev_for_prod_exports = "Forest Product Exports",
    post_rta = "Standard RTA Effect",
    "post_rta:enviro_rta" = "Enviro Provision Effect"
  ),
  
  fixef_sizes = FALSE,
  
  tex = TRUE,
  file = out_file,
  replace = FALSE,
  title = NULL
  # notes = "Note: All regressions include country and year fixed effects. Robust covariance matrices accounting for cross-cluster correlation are used for all specifications."
)

cat("}\n\n\\vspace{2em}\n\n", file = out_file, append = TRUE)


############################################################
# Regressions on tropical "tropical_harvest_ha","tropical_harvest_ton","tropical_harvest_yield","tropical_ag_exp_val", "tropical_ag_exp_val"/"tropical_ag_exp_ton","tropical_for_prod_exports", "for_for_prod_output"

# Regressions on tropical country variables (matched data only)
matched_data$ln_tropical_harvest_ha <- log(matched_data$tropical_harvest_ha + 1)
matched_data$ln_tropical_harvest_ton <- log(matched_data$tropical_harvest_ton + 1)
matched_data$ln_tropical_harvest_yield <- log(matched_data$tropical_harvest_yield + 1)
matched_data$ln_tropical_ag_exp_val <- log(matched_data$tropical_ag_exp_val + 1)
matched_data$ln_tropical_ag_exp_unit_val <- log(
  matched_data$tropical_ag_exp_val / matched_data$tropical_ag_exp_ton + 1)
  matched_data[is.na(matched_data)] <- 0
matched_data$ln_tropical_for_prod_exports <- log(
  matched_data$tropical_for_prod_exports + 1)
  matched_data$ln_tropical_for_prod_output <- log(
    matched_data$tropical_for_prod_output + 1)


tropical_harvest_ha_model <- feols(
  ln_tropical_harvest_ha ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_tropical_harvest_ha <- compute_robust_vcov(tropical_harvest_ha_model, W_matrix, matched_data, unique_ids)

tropical_harvest_ton_model <- feols(
  ln_tropical_harvest_ton ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_tropical_harvest_ton <- compute_robust_vcov(tropical_harvest_ton_model, W_matrix, matched_data, unique_ids)

tropical_harvest_yield_model <- feols(
  ln_tropical_harvest_yield ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_tropical_harvest_yield <- compute_robust_vcov(tropical_harvest_yield_model, W_matrix, matched_data, unique_ids)

tropical_ag_exp_val_model <- feols(
  ln_tropical_ag_exp_val ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_tropical_harvest_ag_exp_val <- compute_robust_vcov(tropical_ag_exp_val_model, W_matrix, matched_data, unique_ids)

tropical_ag_exp_unit_val_model <- feols(
  ln_tropical_ag_exp_unit_val ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_tropical_ag_exp_unit_val <- compute_robust_vcov(tropical_ag_exp_unit_val_model, W_matrix, matched_data, unique_ids)

  tropical_for_prod_exports_model <- feols(
    ln_tropical_for_prod_exports ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_tropical_for_prod_exports <- compute_robust_vcov(tropical_for_prod_exports_model, W_matrix, matched_data, unique_ids)

    tropical_for_prod_outputs_model <- feols(
      ln_tropical_for_prod_output ~ post_rta + post_rta:enviro_rta | id + year,
      data = matched_data,
      cluster = ~id,
      weights = ~weights)

    V_tropical_for_prod_output <- compute_robust_vcov(tropical_for_prod_outputs_model, W_matrix, matched_data, unique_ids)

tropical_for_prod_exports_model <- feols(
  ln_tropical_for_prod_exports ~ post_rta + post_rta:enviro_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_tropical_for_prod_exports <- compute_robust_vcov(tropical_for_prod_exports_model, W_matrix, matched_data, unique_ids)

# Table 4b: Mechanism Analysis - Tropical Countries
cat("\\section*{Table 4B: Mechanism - Tropical Countries Agriculture \\& Forestry}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

etable(
  tropical_harvest_ha_model,
  tropical_harvest_ton_model,
  tropical_harvest_yield_model,
  tropical_ag_exp_val_model,
  tropical_ag_exp_unit_val_model,
  tropical_for_prod_outputs_model,
  tropical_for_prod_exports_model,
  signif.code = c("***"=0.01, "**"=0.05, "*"=0.10),
  vcov = list(
    V_tropical_harvest_ha,
    V_tropical_harvest_ton,
    V_tropical_harvest_yield,
    V_tropical_harvest_ag_exp_val,
    V_tropical_ag_exp_unit_val,
    V_tropical_for_prod_output,
    V_tropical_for_prod_exports
  ),
  dict = c(
    ln_tropical_harvest_ha = "Harvest Area (ha)",
    ln_tropical_harvest_ton = "Harvest (tons)",
    ln_tropical_harvest_yield = "Harvest Yield",
    ln_tropical_ag_exp_val = "Ag Export Value",
    ln_tropical_ag_exp_unit_val = "Ag Export Unit Value",
    tropical_for_prod_output = "Forest Product Output",
    tropical_for_prod_exports = "Forest Product Exports",
    post_rta = "Standard RTA Effect",
    "post_rta:enviro_rta" = "Enviro Provision Effect"
  ),
  
  fixef_sizes = FALSE,
  
  tex = TRUE,
  file = out_file,
  replace = FALSE,
  title = NULL
  # notes = "Note: All regressions include country and year fixed effects. Robust covariance matrices accounting for cross-cluster correlation are used for all specifications."
  )

cat("}\n\n\\vspace{2em}\n\n", file = out_file, append = TRUE)


##########################################################
# Regressions on government expenditure variables for developing and tropical countries
  ##########################################################
  # Matched regression for developing countries
  dev_ag_recur_tot_model <- feols(
    log(dev_ag_recur_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_dev_ag_recur_tot <- compute_robust_vcov(dev_ag_recur_tot_model, W_matrix, matched_data, unique_ids)

  dev_rnd_ag_tot_model <- feols(
    log(dev_rnd_ag_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_dev_rnd_ag_tot <- compute_robust_vcov(dev_rnd_ag_tot_model, W_matrix, matched_data, unique_ids)
    
    # remove negative
    matched_data$dev_ag_cap_tot[matched_data$dev_ag_cap_tot < 0] <- 0
  dev_ag_cap_tot_model <- feols(
    log(dev_ag_cap_tot+1)  ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_dev_ag_cap_tot <- compute_robust_vcov(dev_ag_cap_tot_model, W_matrix, matched_data, unique_ids)


    # remove negative
    matched_data$dev_for_recur_tot[matched_data$dev_for_recur_tot < 0] <- 0
  dev_for_recur_tot_model <- feols(
    log(dev_for_recur_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_dev_for_recur_tot <- compute_robust_vcov(dev_for_recur_tot_model, W_matrix, matched_data, unique_ids)
    

    # remove negative
    matched_data$dev_for_cap_tot[matched_data$dev_for_cap_tot < 0] <- 0
  dev_for_cap_tot_model <- feols(
    log(dev_for_cap_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_dev_for_cap_tot <- compute_robust_vcov(dev_for_cap_tot_model, W_matrix, matched_data, unique_ids)

  # Table 5A: Government Expenditure - Developing Countries
  cat("\\section*{Table 5A: Government Expenditure - Developing Countries}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

  etable(
    dev_ag_recur_tot_model,
    dev_rnd_ag_tot_model,
    dev_ag_cap_tot_model,
    dev_for_recur_tot_model,
    dev_for_cap_tot_model,
    signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    vcov = list(
      V_dev_ag_recur_tot,
      V_dev_rnd_ag_tot,
      V_dev_ag_cap_tot,
      V_dev_for_recur_tot,
      V_dev_for_cap_tot
    ),
    dict = c(
      dev_ag_recur_tot = "Dev Ag Recurrent Expenditure",
      dev_rnd_ag_tot = "Dev R&D Ag Expenditure",
      dev_ag_cap_tot = "Dev Ag Capital Expenditure",
      dev_for_recur_tot = "Dev Forest Recurrent Expenditure",
      dev_for_cap_tot = "Dev Forest Capital Expenditure",
      post_rta = "Standard RTA Effect",
      "post_rta:enviro_rta" = "Enviro Provision Effect"
    ),
    fixef_sizes = FALSE,
    tex = TRUE,
    file = out_file,
    replace = FALSE,
    title = NULL
  )

cat("}\n\n\\vspace{2em}\n\n", file = out_file, append = TRUE)

  ##########################################################
  # Matched regression for tropical countries
  tropical_ag_recur_tot_model <- feols(
    log(tropical_ag_recur_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_tropical_ag_recur_tot <- compute_robust_vcov(tropical_ag_recur_tot_model, W_matrix, matched_data, unique_ids)

  tropical_rnd_ag_tot_model <- feols(
    log(tropical_rnd_ag_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_tropical_rnd_ag_tot <- compute_robust_vcov(tropical_rnd_ag_tot_model, W_matrix, matched_data, unique_ids)

    # remove negative
    matched_data$tropical_ag_cap_tot[matched_data$tropical_ag_cap_tot < 0] <- 0
  tropical_ag_cap_tot_model <- feols(
    log(tropical_ag_cap_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)


    V_tropical_ag_cap_tot <- compute_robust_vcov(tropical_ag_cap_tot_model, W_matrix, matched_data, unique_ids)

    # remove negative
    matched_data$tropical_for_recur_tot[matched_data$tropical_for_recur_tot < 0] <- 0
  tropical_for_recur_tot_model <- feols(
    log(tropical_for_recur_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_tropical_for_recur_tot <- compute_robust_vcov(tropical_for_recur_tot_model, W_matrix, matched_data, unique_ids)

    # remove negative
    matched_data$tropical_for_cap_tot[matched_data$tropical_for_cap_tot < 0] <- 0

  tropical_for_cap_tot_model <- feols(
    log(tropical_for_cap_tot + 1) ~ post_rta + post_rta:enviro_rta | id + year,
    data = matched_data,
    cluster = ~id,
    weights = ~weights)

    V_tropical_for_cap_tot <- compute_robust_vcov(tropical_for_cap_tot_model, W_matrix, matched_data, unique_ids)

  # Table 5B: Government Expenditure - Tropical Countries
  cat("\\section*{Table 5B: Government Expenditure - Tropical Countries}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

  etable(
    tropical_ag_recur_tot_model,
    tropical_rnd_ag_tot_model,
    tropical_ag_cap_tot_model,
    tropical_for_recur_tot_model,
    tropical_for_cap_tot_model,
    signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    vcov = list(
      V_tropical_ag_recur_tot,
      V_tropical_rnd_ag_tot,
      V_tropical_ag_cap_tot,
      V_tropical_for_recur_tot,
      V_tropical_for_cap_tot
    ),
    dict = c(
      tropical_ag_recur_tot = "Tropical Ag Recurrent Expenditure",
      tropical_rnd_ag_tot = "Tropical R&D Ag Expenditure",
      tropical_ag_cap_tot = "Tropical Ag Capital Expenditure",
      tropical_for_recur_tot = "Tropical Forest Recurrent Expenditure",
      tropical_for_cap_tot = "Tropical Forest Capital Expenditure",
      post_rta = "Standard RTA Effect",
      "post_rta:enviro_rta" = "Enviro Provision Effect"
    ),
    
    fixef_sizes = FALSE,
    
    tex = TRUE,
    file = out_file,
    replace = FALSE,
    title = NULL
    # notes = "Note: All regressions include country and year fixed effects. Robust covariance matrices accounting for cross-cluster correlation are used for all specifications."
  )

cat("}\n\n\\vspace{2em}\n\n", file = out_file, append = TRUE)





############################################################
# Matched and non-matched regression of log forest loss and deforestation rate using post_rta, enviro_rta and enviro_enforce_rta 


# Regression with ln_loss (matched)
forest_prov_model <- feols(
  ln_loss ~ post_rta + post_rta:enviro_rta + enviro_enforce_rta |
    id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

  V_forest_prov <- compute_robust_vcov(forest_prov_model, W_matrix, matched_data, unique_ids)

# Regression with ln_loss (not matched)
forest_prov_model_not_matched <- feols(
  ln_loss ~ post_rta + post_rta:enviro_rta + enviro_enforce_rta |
    id + year,
  data = panel_data,
  cluster = ~id)

  V_forest_prov_not_matched <- compute_robust_vcov(forest_prov_model_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)

# Regression with average loss rate (matched)
forest_prov_rate_model <- feols(
  avg_loss_rate ~ post_rta + post_rta:enviro_rta +
    enviro_enforce_rta | id + year,
  data = matched_data,
  cluster = ~id,
  weights = ~weights)

    V_forest_prov_rate <- compute_robust_vcov(forest_prov_rate_model, W_matrix, matched_data, unique_ids)



# Regression with average loss rate (not matched)
forest_prov_rate_model_not_matched <- feols(
  avg_loss_rate ~ post_rta + post_rta:enviro_rta +
    enviro_enforce_rta | id + year,
  data = panel_data,
  cluster = ~id)

    V_forest_prov_rate_not_matched <- compute_robust_vcov(forest_prov_rate_model_not_matched, W_matrix_not_matched, panel_data, unique_ids_panel)


# Table 6: Environmental Enforcement Effects
cat("\\section*{Table 6: Environmental Enforcement Effects}

\\resizebox{\\textwidth}{!}{%
", file = out_file, append = TRUE)

etable(
  forest_prov_model,
  forest_prov_model_not_matched,
  forest_prov_rate_model,
  forest_prov_rate_model_not_matched,
  signif.code = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
  vcov = list(
    V_forest_prov,
    V_forest_prov_not_matched,
    V_forest_prov_rate,
    V_forest_prov_rate_not_matched
  ),
  dict = c(
    ln_loss = "Log Forest Loss",
    avg_loss_rate = "Average Loss Rate",
    post_rta = "Standard RTA Effect",
    "post_rta:enviro_rta" = "Enviro Provision Effect",
    enviro_enforce_rta = "Enviro Enforcement Effect"
  ),
  
  fixef_sizes = FALSE,
  tex = TRUE,
  file = out_file,
  replace = FALSE,
  title = NULL
  # notes = "Note: All regressions include country and year fixed effects. Columns 1 and 3 use robust covariance matrix accounting for cross-cluster correlation (matched sample). Columns 2 and 4 use standard errors clustered at country level (full panel)."
)

cat("}\n\n\\end{document}", file = out_file, append = TRUE)


