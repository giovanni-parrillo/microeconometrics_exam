library(fixest)
library(lubridate)
library(dplyr)

# 1. Load the fixed dataset
# Make sure this is the file you just created with the post-years included
matched_data <- read.csv("matched_data_fixed.csv")

# 2. Robust Date Parsing
# We try multiple formats: "y" (2005), "ymd" (2005-01-01), "mdy" (01/01/2005)
matched_data$entry_date_parsed <- parse_date_time(matched_data$Entry.into.Force, 
                                                  orders = c("y", "ymd", "dmy", "mdy"))

# Extract the numeric year
matched_data$entry_year <- year(matched_data$entry_date_parsed)

# 3. Create Variables
# Create Post variable (1 if year >= entry_year)
matched_data$post_rta <- ifelse(matched_data$year >= matched_data$entry_year, 1, 0)

# Create Log Loss (Outcome)
matched_data$ln_loss <- log(matched_data$loss + 1)

# Handle weights (fill NAs with 0)
matched_data$weights[is.na(matched_data$weights)] <- 0

# 4. Final Diagnostic Check
# You must see numbers in all 4 quadrants (including the top-right!)
print("--- Distribution check ---")
print(table(Enviro = matched_data$enviro_rta, Post = matched_data$post_rta))

# 5. Run Model
model <- feols(ln_loss ~ post_rta + post_rta:enviro_rta | id + year,
               data = matched_data,
               cluster = ~id,
               weights = ~weights)

# 6. View Results
etable(model, 
       signifCode = c("***"=0.01, "**"=0.05, "*"=0.10),
       dict = c(ln_loss = "Log Forest Loss",
                post_rta = "Standard RTA Effect",
                "post_rta:enviro_rta" = "Enviro Provision Effect"))
