# =============================================================================
# Optimized Estimation Pipeline with Robust Covariance Matrix
# =============================================================================

library(fixest)
library(dplyr)

# =============================================================================
# CORE ESTIMATION FUNCTION
# =============================================================================

robust_feols <- function(data,
                         dep_var,
                         indep_vars = "post_rta + post_rta:enviro_rta",
                         fixed_effects = "id + year",
                         cluster_var = "id",
                         rta_cols = 2:202,
                         use_weights = FALSE,
                         display = TRUE,
                         return_latex = FALSE) {
  #' Estimate fixed-effects model with robust covariance matrix

  #' 
  #' @param data Dataset for estimation
  #' @param dep_var Dependent variable name (string)
  #' @param indep_vars Independent variables specification (string)
  #' @param fixed_effects Fixed effects specification (string)
  #' @param cluster_var Variable name for clustering (string)
  #' @param rta_cols Column indices for RTA membership indicators
  #' @param use_weights Whether to use propensity score weights
  #' @param display Whether to display results
  #' @param return_latex Whether to return LaTeX output
  #' @return List with model, robust vcov, and optionally LaTeX output
  
  # --- Step 1: Build weight matrix ---
  unique_ids <- unique(data[[cluster_var]])
  n_clusters <- length(unique_ids)
  
  cat("Building weight matrix for", n_clusters, "clusters...\n")
  
  W <- matrix(0, n_clusters, n_clusters)
  
  for (i in seq_len(n_clusters)) {
    g_data <- data[data[[cluster_var]] == unique_ids[i], ][1, rta_cols]
    sum_g <- sum(g_data, na.rm = TRUE)
    
    if (sum_g > 0) {
      for (k in seq_len(n_clusters)) {
        h_data <- data[data[[cluster_var]] == unique_ids[k], ][1, rta_cols]
        both_ones <- sum((g_data == 1) & (h_data == 1), na.rm = TRUE)
        W[i, k] <- both_ones / sum_g
      }
    }
  }
  
  cat("Weight matrix complete.\n")
  
  # --- Step 2: Estimate model ---
  formula_str <- paste(dep_var, "~", indep_vars, "|", fixed_effects)
  formula_obj <- as.formula(formula_str)
  cluster_formula <- as.formula(paste0("~", cluster_var))
  
  cat("Estimating model:", formula_str, "\n")
  
  if (use_weights && "weights" %in% names(data)) {
    model <- feols(formula_obj, data = data, cluster = cluster_formula, weights = ~weights)
  } else {
    model <- feols(formula_obj, data = data, cluster = cluster_formula)
  }
  
  # --- Step 3: Compute robust covariance matrix ---
  eps <- residuals(model)
  X <- model.matrix(model)
  
  # Observation-level scores
  obs_scores <- sweep(X, 1, eps, '*')
  
  # Aggregate by cluster
  group_indices <- match(data[[cluster_var]], unique_ids)
  S <- rowsum(obs_scores, group = group_indices)
  
  # Sandwich formula: V = (X'X)^{-1} * S'WS * (X'X)^{-1}
  XtX_inv <- solve(crossprod(X))
  V <- XtX_inv %*% crossprod(S, W %*% S) %*% XtX_inv
  
  cat("Robust covariance matrix computed.\n\n")
  
  # --- Step 4: Display results ---
  if (display) {
    cat("========================================\n")
    cat("Results for:", dep_var, "\n")
    cat("========================================\n")
    
    print(etable(model,
                 vcov = V,
                 signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
                 dict = c(post_rta = "Post-RTA",
                          "post_rta:enviro_rta" = "Env. Provision Effect")))
  }
  
  # --- Step 5: Prepare output ---
  result <- list(
    model = model,
    vcov = V,
    W_matrix = W,
    unique_ids = unique_ids,
    formula = formula_str
  )
  
  if (return_latex) {
    result$latex <- etable(model,
                           vcov = V,
                           signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
                           dict = c(post_rta = "Post-RTA",
                                    "post_rta:enviro_rta" = "Env. Provision Effect"),
                           tex = TRUE)
  }
  
  return(result)
}

# =============================================================================
# COMPARISON FUNCTION (MATCHED VS PANEL)
# =============================================================================

compare_samples <- function(matched_data,
                            panel_data,
                            dep_var,
                            indep_vars = "post_rta + post_rta:enviro_rta",
                            fixed_effects = "id + year",
                            rta_cols = 2:202,
                            export_latex = NULL) {
  #' Compare estimation results between matched and full panel samples
  #' 
  #' @param matched_data Matched sample dataset
  #' @param panel_data Full panel dataset
  #' @param dep_var Dependent variable name
  #' @param indep_vars Independent variables specification
  #' @param fixed_effects Fixed effects specification
  #' @param rta_cols Column indices for RTA membership
  #' @param export_latex Filename for LaTeX export (NULL = no export)
  #' @return List with both estimation results
  
  cat("\n============================================\n")
  cat("  MATCHED SAMPLE ESTIMATION\n")
  cat("============================================\n\n")
  
  res_matched <- robust_feols(
    data = matched_data,
    dep_var = dep_var,
    indep_vars = indep_vars,
    fixed_effects = fixed_effects,
    rta_cols = rta_cols,
    use_weights = TRUE,
    display = TRUE
  )
  
  cat("\n============================================\n")
  cat("  FULL PANEL ESTIMATION\n")
  cat("============================================\n\n")
  
  res_panel <- robust_feols(
    data = panel_data,
    dep_var = dep_var,
    indep_vars = indep_vars,
    fixed_effects = fixed_effects,
    rta_cols = rta_cols,
    use_weights = FALSE,
    display = TRUE
  )
  
  # Combined table
  cat("\n============================================\n")
  cat("  COMBINED COMPARISON TABLE\n")
  cat("============================================\n\n")
  
  combined_table <- etable(
    res_matched$model, res_panel$model,
    vcov = list(res_matched$vcov, res_panel$vcov),
    signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    headers = c("Matched", "Full Panel"),
    dict = c(post_rta = "Post-RTA",
             "post_rta:enviro_rta" = "Env. Provision Effect")
  )
  
  print(combined_table)
  
  # Export to LaTeX if requested
  if (!is.null(export_latex)) {
    etable(
      res_matched$model, res_panel$model,
      vcov = list(res_matched$vcov, res_panel$vcov),
      signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
      headers = c("Matched", "Full Panel"),
      dict = c(post_rta = "Post-RTA",
               "post_rta:enviro_rta" = "Env. Provision Effect"),
      title = paste("Effect on", dep_var),
      notes = "Robust SE with cross-cluster correlation in parentheses.",
      tex = TRUE,
      file = export_latex,
      style.tex = style.tex("aer"),
      fitstat = c("n", "r2")
    )
    cat("\nLaTeX exported to:", export_latex, "\n")
  }
  
  return(list(matched = res_matched, panel = res_panel))
}

# =============================================================================
# UTILITY: REUSE WEIGHT MATRIX
# =============================================================================

robust_feols_with_W <- function(data,
                                dep_var,
                                W_matrix,
                                unique_ids,
                                indep_vars = "post_rta + post_rta:enviro_rta",
                                fixed_effects = "id + year",
                                cluster_var = "id",
                                use_weights = FALSE,
                                display = TRUE) {
  #' Estimate model reusing a pre-computed weight matrix
  #' Use this for multiple dependent variables on the same dataset
  #' 
  #' @param data Dataset
  #' @param dep_var Dependent variable name
  #' @param W_matrix Pre-computed weight matrix
  #' @param unique_ids Vector of unique cluster IDs matching W_matrix
  #' @param indep_vars Independent variables
  #' @param fixed_effects Fixed effects
  #' @param cluster_var Cluster variable name
  #' @param use_weights Use propensity weights
  #' @param display Display results
  #' @return List with model and vcov
  
  # Estimate model
  formula_str <- paste(dep_var, "~", indep_vars, "|", fixed_effects)
  formula_obj <- as.formula(formula_str)
  cluster_formula <- as.formula(paste0("~", cluster_var))
  
  if (use_weights && "weights" %in% names(data)) {
    model <- feols(formula_obj, data = data, cluster = cluster_formula, weights = ~weights)
  } else {
    model <- feols(formula_obj, data = data, cluster = cluster_formula)
  }
  
  # Compute robust vcov
  eps <- residuals(model)
  X <- model.matrix(model)
  obs_scores <- sweep(X, 1, eps, '*')
  group_indices <- match(data[[cluster_var]], unique_ids)
  S <- rowsum(obs_scores, group = group_indices)
  XtX_inv <- solve(crossprod(X))
  V <- XtX_inv %*% crossprod(S, W_matrix %*% S) %*% XtX_inv
  
  if (display) {
    cat("\n--- Results for:", dep_var, "---\n")
    print(etable(model, vcov = V,
                 signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
                 dict = c(post_rta = "Post-RTA",
                          "post_rta:enviro_rta" = "Env. Provision Effect")))
  }
  
  return(list(model = model, vcov = V))
}


# =============================================================================
# FASTER COMPARISON
# =============================================================================
compare_samples_with_W <- function(matched_data,
                            panel_data,
                            dep_var,
                            indep_vars = "post_rta + post_rta:enviro_rta",
                            fixed_effects = "id + year",
                            rta_cols = 2:202,
                            export_latex = NULL) {
  #' Compare estimation results between matched and full panel samples
  #' 
  #' @param matched_data Matched sample dataset
  #' @param panel_data Full panel dataset
  #' @param dep_var Dependent variable name
  #' @param indep_vars Independent variables specification
  #' @param fixed_effects Fixed effects specification
  #' @param rta_cols Column indices for RTA membership
  #' @param export_latex Filename for LaTeX export (NULL = no export)
  #' @return List with both estimation results
  
  cat("\n============================================\n")
  cat("  MATCHED SAMPLE ESTIMATION\n")
  cat("============================================\n\n")
  
  res_matched <- robust_feols_with_W(
    data = matched_data,
    dep_var = dep_var,
    W_matrix = res1$W_matrix,
    unique_ids = res1$unique_ids,
    indep_vars = indep_vars,
    fixed_effects = fixed_effects,
    rta_cols = rta_cols,
    use_weights = TRUE,
    display = TRUE
  )
  
  cat("\n============================================\n")
  cat("  FULL PANEL ESTIMATION\n")
  cat("============================================\n\n")
  
  res_panel <- robust_feols_with_W(
    data = panel_data,
    dep_var = dep_var,
    W_matrix = res1$W_matrix,
    unique_ids = res1$unique_ids,
    indep_vars = indep_vars,
    fixed_effects = fixed_effects,
    rta_cols = rta_cols,
    use_weights = FALSE,
    display = TRUE
  )
  
  # Combined table
  cat("\n============================================\n")
  cat("  COMBINED COMPARISON TABLE\n")
  cat("============================================\n\n")
  
  combined_table <- etable(
    res_matched$model, res_panel$model,
    vcov = list(res_matched$vcov, res_panel$vcov),
    signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
    headers = c("Matched", "Full Panel"),
    dict = c(post_rta = "Post-RTA",
             "post_rta:enviro_rta" = "Env. Provision Effect")
  )
  
  print(combined_table)
  
  # Export to LaTeX if requested
  if (!is.null(export_latex)) {
    etable(
      res_matched$model, res_panel$model,
      vcov = list(res_matched$vcov, res_panel$vcov),
      signifCode = c("***" = 0.01, "**" = 0.05, "*" = 0.10),
      headers = c("Matched", "Full Panel"),
      dict = c(post_rta = "Post-RTA",
               "post_rta:enviro_rta" = "Env. Provision Effect"),
      title = paste("Effect on", dep_var),
      notes = "Robust SE with cross-cluster correlation in parentheses.",
      tex = TRUE,
      file = export_latex,
      style.tex = style.tex("aer"),
      fitstat = c("n", "r2")
    )
    cat("\nLaTeX exported to:", export_latex, "\n")
  }
  
  return(list(matched = res_matched, panel = res_panel))
}




# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# # Load data
matched_data <- read.csv("matched_data_fixed.csv")
panel_data <- read.csv("rta_panel.csv")

# Prepare variables
matched_data$post_rta <- ifelse(matched_data$year >= matched_data$Entry.into.Force, 1, 0)
matched_data$ln_loss <- log(matched_data$loss + 1)
matched_data$weights[is.na(matched_data$weights)] <- 0

panel_data$post_rta <- ifelse(panel_data$year >= panel_data$Entry.into.Force, 1, 0)
panel_data$ln_loss <- log(panel_data$loss + 1)

# ----- Compare matched vs panel -----
results <- compare_samples(matched_data, panel_data, "ln_loss", 
                           export_latex = "results.tex")

# ----- Multiple dep vars (reusing weight matrix) -----
# First run computes W
res1 <- robust_feols(matched_data, "ln_loss", use_weights = TRUE)

# Subsequent runs reuse W (much faster)
res2 <- robust_feols_with_W(matched_data, "avg_loss_rate", 
                             res1$W_matrix, res1$unique_ids, use_weights = TRUE)



# Changing dependent and independent variables

# ----- Single estimation -----
result <- robust_feols(matched_data, "ln_loss", use_weights = TRUE)

# ----- Compare matched vs panel -----
results <- compare_samples(matched_data, panel_data, "ln_loss", 
                           export_latex = "results.tex")

# ----- Multiple dep vars (reusing weight matrix) -----
# First run computes W
res1 <- robust_feols(matched_data, "ln_loss", use_weights = TRUE)

# Subsequent runs reuse W (much faster)
res2 <- robust_feols_with_W(matched_data, "avg_loss_rate", 
                            res1$W_matrix, res1$unique_ids, use_weights = TRUE)

ciccio <- compare_samples_with_W(matched_data, panel_data, "ln_loss", 
                                 export_latex = "results2.tex")



