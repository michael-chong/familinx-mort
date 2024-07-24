library(tidyverse)
library(here)

for (subset_name in str_c("subset_", 1:9)) {
  
  #prof <- read_csv(here(paste0("data/intermediate/", subset_name, ".csv")))
  prof <- read_rds(here(paste0("data/intermediate/", subset_name, ".rds")))
  
  # Easy fixes
  prof <- prof %>%
    select(-ends_with("2")) %>% # these extra location columns are redundant
    select(-contains("latitude"), -contains("longitude")) %>% # don't have the time to geocode locations
    select(-is_alive) %>% # for our time period, individuals shouldn't still be alive
    select(-contains("extern")) %>%
    select(-contains("day"), -contains("month")) # too granular for our purposes
  
  # Date stuff
  prof <- prof %>% 
    # Try to impute birth dates with baptism dates
    mutate(birth_year = ifelse(
      is.na(birth_year),
      baptism_year,
      birth_year
    )) %>%
    mutate(birth_date_text = ifelse(
      is.na(birth_date_text),
      baptism_date_text,
      birth_date_text
    )) %>%
    # Try to impute death dates with burial dates
    mutate(death_year = ifelse(
      is.na(death_year),
      burial_year,
      death_year
    )) %>%
    mutate(death_date_text = ifelse(
      is.na(death_date_text),
      burial_date_text,
      death_date_text
    )) %>%
    # Parse year from text fields if not present
    mutate(birth_year = ifelse(
      is.na(birth_year),
      as.numeric(str_extract(birth_date_text, "[0-9][0-9][0-9][0-9]")), 
      birth_year)) %>%
    mutate(death_year = ifelse(
      is.na(death_year),
      as.numeric(str_extract(death_date_text, "[0-9][0-9][0-9][0-9]")),
      death_year))
  
  # Location stuff
  prof <- prof %>%
    # Replace NAs with blank strings
    mutate(across(contains("location"), ~str_replace_na(.x, ""))) %>%
    # Combine city, state, and country to one string (sometimes country is in the wrong column)
    mutate(birth_place = str_c(birth_location_city, birth_location_state, birth_location_country, sep = ", ")) %>%
    mutate(baptism_place = str_c(baptism_location_city, baptism_location_state, baptism_location_country, sep = ", ")) %>%
    mutate(death_place = str_c(death_location_city, death_location_state, death_location_country, sep = ", ")) %>%
    mutate(burial_place = str_c(burial_location_city, burial_location_state, burial_location_country, sep = ", ")) %>%
    # Impute places with the free text field
    mutate(birth_place = ifelse(
      is.na(birth_place),
      birth_location_place_name,
      birth_place)) %>%
    mutate(baptism_place = ifelse(
      is.na(baptism_place),
      baptism_location_place_name,
      baptism_place)) %>%
    mutate(death_place = ifelse(
      is.na(death_place),
      death_location_place_name,
      death_place)) %>%
    mutate(burial_place = ifelse(
      is.na(burial_place),
      burial_location_place_name,
      burial_place)) %>%
    # Impute birth place with baptism place
    mutate(birth_place = ifelse(
      is.na(birth_place),
      baptism_place,
      birth_place)) %>%
    # Impute death place with burial place
    mutate(death_place = ifelse(
      is.na(death_place),
      burial_place,
      death_place))
  
  prof <- prof %>% select(profileid, gender, birth_year, death_year, birth_place, death_place)
  
  write_csv(prof, here(paste0("data/intermediate/", subset_name, "_narrowed.csv")))
}

