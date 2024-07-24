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
  mutate(birth_year = start_year, death_year = end_year) |>
  as.data.table()

trajectories <- rbindlist(list(parents, non_parents), use.names = TRUE)

certain_match <- future_map_dfr(
  .x = set_names(year_range),
  .id = "year",
  .f = function(y) {
    trajectories |>
      filter(start_year <= y, end_year >= y) |>
      mutate(age = y - birth_year) |>
      filter(start_place == end_place, start_place %in% countries) |> 
      group_by(start_place, age, gender) |>
      summarize(n = n(), .groups = "drop") |>
      as.data.table()
  }) |>
  rename(country = start_place) |>
  mutate(type = "full match")

certain_na <- map_dfr(
  .x = set_names(year_range),
  .id = "year",
  .f = function(y) {
    trajectories |>
      filter(start_year <= y, end_year >= y) |>
      mutate(age = y - birth_year) |>
      filter(
        start_place %in% countries | end_place %in% countries,
        is.na(start_place) | is.na(end_place),
      ) |>
      mutate(country = ifelse(
        is.na(start_place),
        end_place,
        start_place
      )) |>
      group_by(country, age, gender) |>
      summarize(n = n(), .groups = "drop") |>
      as.data.table()
  }) |> 
  mutate(type = "NA endpoint")

origin <- map_dfr(
  .x = set_names(year_range),
  .id = "year",
  .f = function(y) {
    trajectories |>
      filter(start_year <= y, end_year >= y) |>
      mutate(age = y - birth_year) |>
      filter(
        start_place != end_place,
        start_place %in% countries
      ) |>
      rename(country = start_place) |> 
      group_by(country, age, gender) |>
      summarize(n = n(), .groups = "drop") |>
      as_tibble() 
  }) |> 
  mutate(type = "emigrant")

destination <- map_dfr(
  .x = set_names(year_range),
  .id = "year",
  .f = function(y) {
    trajectories |>
      filter(start_year <= y, end_year >= y) |>
      mutate(age = y - birth_year) |>
      filter(
        start_place != end_place,
        end_place %in% countries
      ) |>
      rename(country = end_place) |>
      group_by(country, age, gender) |>
      summarize(n = n(), .groups = "drop") |>
      as.data.table() 
  }) |> 
  mutate(type = "immigrant")

exposures <- rbindlist(
  list(certain_match, certain_na, origin, destination)
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

grouped_exposures <- exposures |>
  mutate(
    age = cut(age, breaks = age_breaks, labels = head(age_breaks, -1), include.lowest = TRUE, right = FALSE),
    year = cut(as.numeric(year), breaks = year_breaks, labels = head(year_breaks, -1), include.lowest = FALSE, right = FALSE)
  ) |>
  mutate(type = fct_relevel(type, c("full match", "NA endpoint", "immigrant", "emigrant"))) |>
  group_by(country, age, year, gender, type) |>
  dplyr::summarize(n = sum(n), .groups = "drop")

write_rds(exposures, here("data/clean/migrant_exposures_granular_new.rds"))
write_rds(grouped_exposures, here("data/clean/migrant_exposures_grouped_new.rds"))
