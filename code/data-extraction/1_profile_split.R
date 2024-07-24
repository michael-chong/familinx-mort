# Split the full data into more manageable portions to process 
library(tidyverse)
library(here)

all_profiles <- data.table::fread(
  here("data/familinx-raw/profiles-anon.txt"),
  na.strings = c("*", ""))

total <- nrow(all_profiles)
subset_size <- ceiling(total/10)

cuts <- c(1, (1:9)*subset_size, total)

# Divide into 9 subsets
for (i in 1:1) {
  profile_subset <- all_profiles[cuts[i]:cuts[i+1], ]
  
  write_rds(profile_subset, str_c("data/intermediate/subset_", i, ".rds"))
  #write_csv(profile_subset, str_c("data/intermediate/profile_subset", i, ".csv"))

}
 