library(fixest)
library(ggplot2)

df <- read.csv("country_panel.csv")

country_model <- feols(log(loss)~
                lead_LR+lead_3+lead_2+rta+lag_1+lag_2+lag_3+lag_LR 
                + enviro_lead_LR+enviro_lead_3+enviro_lead_2+enviro_rta+enviro_lag_1  +enviro_lag_2+enviro_lag_3+enviro_lag_LR
                  | ISO3 + year,
                data = df,
                cluster = ~ISO3)

summary(country_model)

# Extract coefficients and standard errors
country_coefs <- coef(country_model)
country_se <- se(country_model)
country_confint <- confint(country_model)

country_results <- data.frame(
  term = names(country_coefs),
  estimate = country_coefs,
  std_error = country_se,
  conf_low = country_confint[, 1],
  conf_high = country_confint[, 2]
)

# Compute sums of coefficients and standard errors for grouped variables
grouped_vars <- list(
  lead_LR_with = c("lead_LR", "enviro_lead_LR"),
  lead_3_with = c("lead_3", "enviro_lead_3"),
  lead_2_with = c("lead_2", "enviro_lead_2"),
  rta_with = c("rta", "enviro_rta"),
  lag_1_with = c("lag_1", "enviro_lag_1"),
  lag_2_with = c("lag_2", "enviro_lag_2"),
  lag_3_with = c("lag_3", "enviro_lag_3"),
  lag_LR_with = c("lag_LR", "enviro_lag_LR")
)

grouped_results <- data.frame(
  group = character(),
  sum_estimate = numeric(),
  sum_std_error = numeric(),
  stringsAsFactors = FALSE
)

for (group_name in names(grouped_vars)) {
  vars <- grouped_vars[[group_name]]
  sum_coef <- sum(country_coefs[vars], na.rm = TRUE)
  sum_se <- sqrt(sum(country_se[vars]^2, na.rm = TRUE))
  
  grouped_results <- rbind(grouped_results, data.frame(
    term = group_name,
    estimate = sum_coef,
    std_error = sum_se,
    stringsAsFactors = FALSE
  ))
}

# Take the first eight rows of country_results and add them to grouped_results
final_results <- rbind(country_results[1:8, 1:3], grouped_results)
#Add a column that takes the value "treatment" for the last 8 rows and "control" for the first 8 rows
final_results$group <- c(rep("control", 8), rep("treatment", 8))
#Create a column that takes -4 if term starts with "lead_LR", -3 if term starts with "lead_3", -2 if term starts with "lead_2", 0 if term starts with "rta", 1 if term starts with "lag_1", 2 if term starts with "lag_2", 3 if term starts with "lag_3", and 4 if term starts with "lag_LR"
final_results$time <- sapply(final_results$term, function(x) {
  if (grepl("^lead_LR", x)) {
    return(-4)
  } else if (grepl("^lead_3", x)) {
    return(-3)
  } else if (grepl("^lead_2", x)) {
    return(-2)
  } else if (grepl("^rta", x)) {
    return(0)
  } else if (grepl("^lag_1", x)) {
    return(1)
  } else if (grepl("^lag_2", x)) {
    return(2)
  } else if (grepl("^lag_3", x)) {
    return(3)
  } else if (grepl("^lag_LR", x)) {
    return(4)
  }
})

#Add two rows with time =-1 and estimate = 0 and std_error = 0 for both groups
final_results <- rbind(final_results, data.frame(
  term = c("before_control", "before_treatment"),
  estimate = c(0, 0),
  std_error = c(0, 0),
  group = c("control", "treatment"),
  time = c(-1, -1)
))

ggplot(final_results, aes(x = time, y = estimate, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = estimate-1.96*std_error, ymax = estimate+1.96*std_error), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red") +
  labs(title = "Total Effects",
       x = "Event Time",
       y = "Log Forest Loss",
       color = "Group") +
  scale_color_manual(values = c("treatment" = "red", "control" = "black")) +
  scale_x_continuous(
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    labels = c("≤-4", "-3", "-2", "-1", "0", "1", "2", "3", "≥4"),
    limits = c(-4, 4)) +
  theme_minimal()


  # Below a computation of cumulative effects (replicating Table 8)

# lincomb <- function(model, terms){
#   b <- coef(model)
#   V <- vcov(model)

#   L <- setNames(rep(0, length(b)), names(b))
#   L[terms] <- 1

#   est <- sum(L * b)
#   se  <- sqrt(as.numeric(t(L) %*% V %*% L))

#   c(estimate = est, se = se)
# }


# terms_no   <- c("rta", "lag_1", "lag_2", "lag_3")
# terms_env  <- c("enviro_rta", "enviro_lag_1", "enviro_lag_2", "enviro_lag_3")
# terms_with <- c(terms_no, terms_env)

# # Log forest loss
# no_loss   <- lincomb(country_model, terms_no)
# with_loss <- lincomb(country_model, terms_with)
# rel_loss  <- lincomb(country_model, terms_env)

# # # Deforestation rate
# # no_rate   <- lincomb(m_rate, terms_no)
# # with_rate <- lincomb(m_rate, terms_with)
# # rel_rate  <- lincomb(m_rate, terms_env)


# # no_loss
# # with_loss
# # rel_loss
