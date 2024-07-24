library(tidyverse)
library(here)
library(data.table)
library(dtplyr)

profiles <- data.table::fread(here("data/intermediate/profiles_recombined_new.csv"))

links <- data.table::fread(here("data/familinx-raw/relations-anon.txt"))

parent_ids <- unique(links$parent)

profiles_with_birth_info <- profiles |> 
  filter(!is.na(birth_place), !is.na(birth_year)) |> 
  pivot_longer(cols = c("parent1", "parent2"), names_to = "parent", values_to = "parentid") |>
  as.data.table()

# NOTE: we might want to allow the birth year to be missing (since it still might have useful place information)

# parent trajectories ----------
parents <- profiles |>
  filter(
    profileid %in% parent_ids,
    !is.na(birth_year),
    !is.na(death_year)
  ) |>
  inner_join( # this removes parents where the child does not have a birthplace recorded
    profiles_with_birth_info, 
    by = c("profileid" = "parentid"),
    suffix = c("", "_child")
  )

countries <- c("USA", "FIN", "NOR", "DNK", "SWE", "FRA", "BEL", "NLD", "GBR", "CHE")
year_range <- 1800:1900

relevant_parents <- parents |> 
  group_by(profileid) |>
  filter(if_any(contains("place"), ~(.x %in% countries))) |>
  filter(
    death_year >= min(year_range), 
    birth_year <= max(year_range)
  ) |>
  ungroup() |>
  as.data.table() 

parent_traj <- relevant_parents |>
  filter(
    birth_year_child >= birth_year, 
    birth_year_child <= death_year
  ) |>
  group_by(profileid, birth_year, death_year, gender) |>
  arrange(birth_year_child) |> 
  summarize(
    start_year = c(birth_year[1], birth_year_child),
    end_year = c(birth_year_child, death_year[1]),
    start_place = c(birth_place[1], birth_place_child),
    end_place = c(birth_place_child, death_place[1]),
    death_indicator = c(rep(0, length(birth_place_child)), 1)
  ) |>
  as.data.table()

fwrite(parent_traj, file = here("data/intermediate/parent_traj_new.csv"))

# non-parent trajectories ------------
non_parents <- profiles |>
  filter(
    !is.na(birth_year),
    !is.na(death_year)
  ) |>
  anti_join(distinct(parent_traj, profileid), by = c("profileid")) |>
  filter(birth_place %in% countries | death_place %in% countries) |>
  filter(death_year >= min(year_range), birth_year <= max(year_range)) |>
  as.data.table()

non_parent_traj <- non_parents |>
  transmute(
    profileid = profileid, 
    gender = gender,
    start_year = birth_year, 
    end_year = death_year,
    start_place = birth_place, 
    end_place = death_place,
    death_indicator = 1
  ) |>
  as.data.table()

fwrite(non_parent_traj, file = here("data/intermediate/non_parent_traj_new.csv"))

