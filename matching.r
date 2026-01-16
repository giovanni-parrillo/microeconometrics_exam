library(MatchIt)
library(dplyr)

# 1. Load Data
# df = One row per RTA (Cross-section)
# panel_df = One row per RTA-Year (Full time series 2001-2018)
df <- read.csv('rta_propensity.csv')
panel_df <- read.csv('rta_panel.csv')

# 2. Match on the CROSS-SECTION (The Groups)
# We match Agreement to Agreement based on their static characteristics
matchit_model <- matchit(enviro_rta ~ predicted_prob, 
                         data = df,       # <--- CHANGE: Use 'df', not 'panel_df'
                         method = "nearest", 
                         replace = TRUE)

summary(matchit_model)

# 3. Extract the Matched Cross-Section
# This creates a list of the ~44 RTAs used in the analysis (36 Treated + Matches)
# and their correct weights.
matched_groups <- match.data(matchit_model)

# 4. Merge Weights into the Full Panel
# We use inner_join to:
#   a) Filter the panel to keep ONLY the matched RTAs
#   b) Bring the 'weights' column into the panel for the regression
final_dataset <- panel_df %>%
  inner_join(matched_groups[, c("id", "weights")], by = "id")

# 5. Verify and Save
# Check that we have post-years now
print(table(final_dataset$enviro_rta, final_dataset$year > final_dataset$Entry.into.Force))

write.csv(final_dataset, "matched_data_fixed.csv", row.names = FALSE)
