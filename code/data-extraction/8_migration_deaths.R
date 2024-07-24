library(tidyverse)
library(here)
library(data.table)
library(dtplyr)
library(furrr)

plan(multisession, workers = 4)

countries <- c("USA", "FIN", "NOR", "DNK", "SWE", "FRA", "BEL", "NLD", "GBR", "CHE")
year_range <- 1800:1900

parents <- fread(here("data/intermediate/parent_traj_new.csv"), na.strings = "")

non_parents <- fread(here("data/intermediate/non_parent_traj_new.csv"), na.strings = "") |>
  mutate(birth_year = start_year, death_year = end_year)

trajectories <- rbindlist(list(parents, non_parents), use.names = TRUE)

certain <- future_map_dfr(
  .x = set_names(year_range),
  .id = "year",
  .f = function(y) {
    trajectories |>
      filter(death_indicator == 1, end_year == y) |>
      mutate(age = y - birth_year) |>
      filter(end_place %in% countries) |> 
      group_by(end_place, age, gender) |>
      summarize(n = n(), .groups = "drop")
  }) |>
  rename(country = end_place) |>
  mutate(type = "death place specified")
  
imputed <- future_map_dfr(
  .x = set_names(year_range),
  .id = "year",
  .f = function(y) {
    trajectories |>
      filter(death_indicator == 1, end_year == y) |>
      mutate(age = y - birth_year) |>
      filter(is.na(end_place), start_place %in% countries) |> 
      group_by(start_place, age, gender) |>
      summarize(n = n(), .groups = "drop")
  }) |>
  rename(country = start_place) |>
  mutate(type = "imputed")

deaths <- rbindlist(
  list(certain, imputed)
) |> 
  filter(age < 110, !is.na(gender)) |>
  complete(year, country, age, gender, type, fill = list(n = 0))

# Year intervals
year_min <- 1800
year_max <- 1900
year_interval_length <- 5
year_breaks <- seq(year_min, year_max, by = year_interval_length)
n_yb <- length(year_breaks)

# Age groups 
age_min <- 5
age_max <- 90
age_interval_length <- 5
age_breaks <- c(0, 1, seq(age_min, age_max, age_interval_length), Inf)
n_ab <- length(age_breaks)

grouped_deaths <- deaths |>
  mutate(
    age = cut(age, breaks = age_breaks, labels = head(age_breaks, -1), include.lowest = TRUE, right = FALSE),
    year = cut(as.numeric(year), breaks = year_breaks, head(year_breaks, -1), include.lowest = FALSE, right = FALSE)
  ) |>
#  mutate(type = fct_relevel(type, c("full match", "NA endpoint", "immigrant", "emigrant"))) |>
  group_by(country, age, year, gender, type) |>
  dplyr::summarize(n = sum(n))

write_rds(deaths, here("data/clean/migrant_deaths_granular_new.rds"))
write_rds(grouped_deaths, here("data/clean/migrant_deaths_grouped_new.rds"))
