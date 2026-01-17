library(fixest)
library(lubridate)
library(dplyr)

# Load the fixed dataset
# Make sure this is the file you just created with the post-years included
matched_data <- read.csv("matched_data_fixed.csv")
panel_data <- read.csv("rta_panel_full.csv")

# Create Variables
# Create Post variable (1 if year >= entry_year)
matched_data$post_rta <- ifelse(matched_data$year >= matched_data$Entry.into.Force, 1, 0)
panel_data$post_rta <- ifelse(panel_data$year >= panel_data$Entry.into.Force, 1, 0)
# Create Log Loss (Outcome)
matched_data$ln_loss <- log(matched_data$loss+1)
panel_data$ln_loss <- log(panel_data$loss+1) #because in panel_data there are losses=0

# Final Diagnostic Check
# You must see numbers in all 4 quadrants (including the top-right!)
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

W_matrix <- matrix(0, length(unique_ids), length(unique_ids))

for (i in 1:length(unique_ids)){
  for (k in 1:length(unique_ids)){
    W_matrix[i,k] <- w_gen(matched_data, unique_ids[i],unique_ids[k])
  }
}


unique_ids_panel <- unique(panel_data$id)

W_matrix_panel <- matrix(0, length(unique_ids_panel), length(unique_ids_panel))

for (i in 1:length(unique_ids_panel)){
  for (k in 1:length(unique_ids_panel)){
    W_matrix_panel[i,k] <- w_gen(panel_data, unique_ids_panel[i],unique_ids_panel[k])
  }
}

# Run Model to get residuals
model <- feols(ln_loss ~ post_rta + post_rta:enviro_rta | id + year,
               data = matched_data,
               cluster = ~id,
               weights = ~weights)

model_not_matched <- feols(ln_loss ~ post_rta + post_rta:enviro_rta | id + year,
                           data = panel_data,
                           cluster = ~id)

eps <- residuals(model)
X <- model.matrix(model)


obs_scores <- sweep(X,1,eps,'*')

group_indices <- match(matched_data$id, unique_ids)
S_matrix <- rowsum(obs_scores, group = group_indices)

V_hat <- t(S_matrix) %*% W_matrix %*% S_matrix

V <- solve(t(X)%*%X) %*% V_hat %*% solve(t(X)%*%X)


eps_not_matched <- residuals(model_not_matched)
X_not_matched <- model.matrix(model_not_matched)


obs_scores_panel <- sweep(X_not_matched,1,eps_not_matched,'*')

group_indices_panel <- match(panel_data$id, unique_ids_panel)
S_matrix_panel <- rowsum(obs_scores_panel, group = group_indices_panel)

V_hat_panel <- t(S_matrix_panel) %*% W_matrix_panel %*% S_matrix_panel

V_panel <- solve(t(X_not_matched)%*%X_not_matched) %*% V_hat_panel %*% solve(t(X_not_matched)%*%X_not_matched)


#View Results
etable(model, 
       signifCode = c("***"=0.01, "**"=0.05, "*"=0.10),
       vcov = V,
       dict = c(ln_loss = "Log Forest Loss",
                post_rta = "Standard RTA Effect",
                "post_rta:enviro_rta" = "Enviro Provision Effect"))

etable(model_not_matched, 
       signifCode = c("***"=0.01, "**"=0.05, "*"=0.10),
       vcov = V_panel,
       dict = c(ln_loss = "Log Forest Loss",
                post_rta = "Standard RTA Effect",
                "post_rta:enviro_rta" = "Enviro Provision Effect"))


#Other meaningful regressions

rate <- feols(avg_loss_rate ~ post_rta + post_rta:enviro_rta | id + year,
                           data = matched_data,
                           cluster = ~id,
                           weights = ~weights)

rate_not_matched <- feols(avg_loss_rate ~ post_rta + post_rta:enviro_rta | id + year,
                            data = panel_data,
                            cluster = ~id)

etable(rate, 
       signifCode = c("***"=0.01, "**"=0.05, "*"=0.10),
       vcov = V,
       dict = c(ln_loss = "Log Forest Loss",
                post_rta = "Standard RTA Effect",
                "post_rta:enviro_rta" = "Enviro Provision Effect"))

etable(rate_not_matched, 
       vcov = V_panel,
       signifCode = c("***"=0.01, "**"=0.05, "*"=0.10),
       dict = c(ln_loss = "Log Forest Loss",
                post_rta = "Standard RTA Effect",
                "post_rta:enviro_rta" = "Enviro Provision Effect"))