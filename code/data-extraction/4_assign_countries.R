library(tidyverse)
library(furrr)
library(here)

plan(multisession, workers = 4)

assign_country <- function(prof, match_expr, country_code) {
  prof %>%
    mutate(birth_place = ifelse(str_detect(birth_place, match_expr), country_code, birth_place)) %>%
    mutate(death_place = ifelse(str_detect(death_place, match_expr), country_code, death_place))
}

country_codes <- read_rds("data/country_codes.rds")

future_walk(
  .x = str_c("subset_", 1:9),
  .f = function(x) {
    profiles <- read_csv(here(paste0("data/intermediate/", x, "_narrowed.csv"))) %>%
      mutate(across(c("birth_place", "death_place"), str_to_upper))
    
    for (country in names(country_codes)) {
      
      if (country != "GBR") {
        profiles <- assign_country(
          profiles,
          match_expr = str_c(country_codes[[country]], collapse = "|"),
          country
        )
      } else {
        
        GBR_match_expr <- str_c(country_codes[[country]], collapse = "|")
        ireland_expr <- c("IRELAND", "DUBLIN", "CORK", "BELFAST", "DERRY", "MUNSTER", "LEINSTER", "ULSTER") |>
          str_c(collapse = "|")
          
        profiles <- profiles |>
          mutate(
            birth_place = ifelse(
              str_detect(birth_place, GBR_match_expr) & !str_detect(birth_place, ireland_expr),
              "GBR",
              birth_place
            ),
            death_place = ifelse(
              str_detect(death_place, GBR_match_expr) & !str_detect(death_place, ireland_expr),
              "GBR",
              death_place
            )
          )
      }

    }
    
    write_csv(profiles, here(paste0("data/intermediate/", x, "_assigned_new.csv")))
  }
)
