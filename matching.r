library(MatchIt)
library(dplyr)

df <- read.csv('rta_propensity.csv')
panel_df <- read.csv('rta_panel_full.csv')

# We match Agreement to Agreement based on their static characteristics
matchit_model <- matchit(enviro_rta ~ predicted_prob, 
                         data = df,
                         method = "nearest", 
                         replace = TRUE)

summary(matchit_model)

# Extract the Matched Cross-Section
matched_groups <- match.data(matchit_model)

# Merge Weights into the Full Panel

final_dataset <- panel_df %>%
  inner_join(matched_groups[, c("id", "weights")], by = "id")

print(table(final_dataset$enviro_rta, final_dataset$year > final_dataset$Entry.into.Force))

write.csv(final_dataset, "matched_data_fixed.csv", row.names = FALSE)
