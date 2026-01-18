library(dplyr)
library(ggplot2)

matched_data <- read.csv("matched_data_fixed.csv")

# Create Variables
# Create Post variable (1 if year >= entry_year)
matched_data$post_rta <- ifelse(matched_data$year >= matched_data$Entry.into.Force, 1, 0)
# Create Log Loss (Outcome)
matched_data$ln_loss <- log(matched_data$loss+1) #because in matched_data there are losses=0
#Create a variable for treatment group
matched_data <- matched_data %>%    
  group_by(id) %>%
  mutate(group = ifelse(any(enviro_rta == 1), "treat", "control")) %>%
  ungroup()

# Create time_to_event variable
matched_data <- matched_data |>
  mutate(event_time = year - Entry.into.Force) |>
  mutate(event_time = ifelse(event_time > 4, 4, ifelse(event_time < -4, -4, event_time))) |>
  mutate(event_time = as.integer(event_time))

matched_data <- dummy_cols(matched_data, select_columns = "event_time")

#Compute average ln_loss by event_time and enviro_rta2, and confidence intervals

ci_ln_loss <- matched_data %>%
  group_by(event_time, group) %>%
  summarize(avg_ln_loss = mean(ln_loss, na.rm = TRUE), 
            sd_ln_loss = sd(ln_loss, na.rm = TRUE),
            n = n(),
            se_ln_loss = sd_ln_loss / sqrt(n),
            ci_lower = avg_ln_loss - qt(0.975, df = n - 1) * se_ln_loss,
            ci_upper = avg_ln_loss + qt(0.975, df = n - 1) * se_ln_loss)


#Plot average ln_loss by event_time and enviro_rta2 with confidence intervals
ggplot(ci_ln_loss, aes(x = event_time, y = avg_ln_loss, color = group)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(title = "Event Study: Impact of RTA on Log Loss",
       x = "Event Time (Years since RTA Entry into Force)",
       y = "Coefficient on Log Loss",
       color = "Group") +
  theme_minimal()

#Now I want to compute average treatment effects by event time
library(fixest)
event_study_model <- feols(ln_loss ~ i(event_time, group, ref = -1) | id + year, data = matched_data)
summary(event_study_model)

#I want to extract the coefficients and confidence intervals for plotting
event_study_coefs <- coef(event_study_model)
event_study_confint <- confint(event_study_model)

# Create a data frame for plotting
# Split row names by "::" to extract components
name_parts <- strsplit(names(event_study_coefs), ":")
event_study_df <- data.frame(
  part1 = sapply(name_parts, function(x) x[1]),
  event_time = as.integer(sapply(name_parts, function(x) ifelse(length(x) >= 3, x[3], NA))),
  part3 = sapply(name_parts, function(x) ifelse(length(x) >= 4, x[4], NA)),
  part4 = sapply(name_parts, function(x) ifelse(length(x) >= 6, x[6], NA)),
  coef = event_study_coefs,
  lower = event_study_confint[, 1],
  upper = event_study_confint[, 2]
)

# Add two reference rows for event_time = -1
ref_row_1 <- data.frame(
  part1 = "event_time",
  event_time = -1,
  part3 = "group",
  part4 = "treat",
  coef = 0,
  lower = 0,
  upper = 0
)
ref_row_2 <- data.frame(
  part1 = "event_time",
  event_time = -1,
  part3 = "group",
  part4 = "control",
  coef = 0,
  lower = 0,
  upper = 0
)
event_study_df <- rbind(ref_row_1, ref_row_2, event_study_df)

# Plotting the event study results
ggplot(event_study_df, aes(x = event_time, y = coef, color = part4)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "red") +
  labs(title = "Event Study: Impact of RTA on Log Loss",
       x = "Event Time (Years since RTA Entry into Force)",
       y = "Coefficient on Log Loss",
       color = "Group") +
  theme_minimal()
