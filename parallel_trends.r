panel_data <- read.csv("rta_panel.csv")

# Create Variables
# Create Post variable (1 if year >= entry_year)
panel_data$post_rta <- ifelse(panel_data$year >= panel_data$Entry.into.Force, 1, 0)
# Create Log Loss (Outcome)
panel_data$ln_loss <- log(panel_data$loss+1) #because in panel_data there are losses=0
#Create a variable for treatment group
panel_data <- panel_data %>%    
  group_by(id) %>%
  mutate(enviro_rta2 = ifelse(any(enviro_rta == 1), "treat", "control")) %>%
  ungroup()
# Final Diagnostic Check
# You must see numbers in all 4 quadrants (including the top-right!)
print("--- Distribution check ---")
print(table(Enviro = panel_data$enviro_rta, Post = panel_data$post_rta))

# Create a new variable that takes value 0 in the first year for which rta==1,
# it takes value 1 in the second year for which rta==1, and so on.
# It should also take negative values for the years before the rta.
panel_data <- panel_data |>
  mutate(event_time = year - Entry.into.Force) |>
  mutate(event_time = ifelse(event_time > 4, 4, ifelse(event_time < -4, -4, event_time)))

panel_data <- dummy_cols(panel_data, select_columns = "event_time")

#Compute average ln_loss by event_time and enviro_rta2, and confidence intervals
library(dplyr)
ci_ln_loss <- panel_data %>%
  group_by(event_time, enviro_rta2) %>%
  summarize(avg_ln_loss = mean(ln_loss, na.rm = TRUE), 
            sd_ln_loss = sd(ln_loss, na.rm = TRUE),
            n = n(),
            se_ln_loss = sd_ln_loss / sqrt(n),
            ci_lower = avg_ln_loss - qt(0.975, df = n - 1) * se_ln_loss,
            ci_upper = avg_ln_loss + qt(0.975, df = n - 1) * se_ln_loss)
print(ci_ln_loss)


#Plot average ln_loss by event_time and enviro_rta2 with confidence intervals
library(ggplot2)
ggplot(ci_ln_loss, aes(x = event_time, y = avg_ln_loss, color = enviro_rta2)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = enviro_rta2), alpha = 0.2, color = NA) +
  labs(title = "Average Log Loss by Event Time and Treatment Group",
       x = "Event Time (Years since RTA Entry into Force)",
       y = "Average Log Loss") +
  theme_minimal()


#Now I want to compute average treatment effects by event time
library(dplyr)
library(fixest)
event_study_model <- feols(ln_loss ~ i(event_time, enviro_rta2, ref = -1) | id2 + year, data = panel_data)
summary(event_study_model)

#I want to extract the coefficients and confidence intervals for plotting
event_study_coefs <- coef(event_study_model)
event_study_confint <- confint(event_study_model)

# Create a data frame for plotting
# Split row names by "::" to extract components
name_parts <- strsplit(names(event_study_coefs), ":")
event_study_df <- data.frame(
  part1 = sapply(name_parts, function(x) x[1]),
  part2 = sapply(name_parts, function(x) ifelse(length(x) >= 2, x[2], NA)),
  part3 = sapply(name_parts, function(x) ifelse(length(x) >= 3, x[3], NA)),
  event_time = as.numeric(gsub("i\\(event_time, enviro_rta2, ref = -1\\)\\[|:.*", "", names(event_study_coefs))),
  coef = event_study_coefs,
  lower = event_study_confint[, 1],
  upper = event_study_confint[, 2]
)



# Plotting the event study results
library(ggplot2)
ggplot(event_study_df, aes(x = event_time, y = coef)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed",
                color = "red") +
    labs(title = "Event Study: Impact of RTA on Log Loss",
         x = "Event Time (Years since RTA Entry into Force)",
         y = "Coefficient on Log Loss") +
    theme_minimal()
