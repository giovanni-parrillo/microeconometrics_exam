library(glmnet)
library(dplyr)
library(fastDummies) # Install if needed: install.packages("fastDummies")
library(fastDummies)

df <- read.csv('rta_data.csv')
df_panel <- read.csv('rta_panel_full.csv')

# Filter to keep only RTAs present in the panel (Consistency check)
df <- df %>% filter(id %in% df_panel$id)

# Ensure 'agree_dev' is a factor and create dummy variables for it
df <- dummy_cols(df, select_columns = "agree_dev", 
                 remove_first_dummy = TRUE, 
                 remove_selected_columns = TRUE)

# Define Predictors (X) and Target (Y)
# We strictly exclude IDs, Outcomes, and Future Characteristics
target_var <- "enviro_rta"

exclude_cols <- c(
  "id", "Agreement", "parties", # IDs / Metadata
  "enviro_rta", "enviro_enforce_rta", # Target & Treatment characteristics
  "defor", "biodiv", "defor_enforce", "biodiv_enforce" # Outcomes (Bad controls)
)

x_data <- df %>% 
  select(where(is.numeric)) %>%
  select(-any_of(exclude_cols)) 

# Convert to Matrix for GLMNET
x_matrix <- data.matrix(x_data)
y_vector <- df[[target_var]]

# Run LASSO
# We use Cross-Validation to find the optimal Lambda
set.seed(123) # For reproducibility
cv_lasso <- cv.glmnet(x_matrix, y_vector, alpha = 1, family = "binomial")

cv_lasso$lambda.min  # Optimal lambda value

# View the coefficients chosen by the model
print(coef(cv_lasso, s = "lambda.min"))

# Predict Propensity Scores
# type="response" gives us the Probability (0 to 1)
probabilities <- predict(cv_lasso, newx = x_matrix, s = "lambda.min", type = "response")

# Add back to dataframe
df$predicted_prob <- as.vector(probabilities)

# Save for the Matching Step
write.csv(df, "rta_propensity.csv", row.names = FALSE)

# If these plots are completely separated, matching will fail.
library(ggplot2)

plot(cv_lasso)

ggplot(df, aes(x = predicted_prob, fill = as.factor(enviro_rta))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Density", fill = "Enviro RTA")
