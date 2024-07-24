library(tidyverse)
library(here)
library(data.table)
library(dtplyr)
library(furrr)

plan(multisession, workers = 4)

links <- data.table::fread(here("data/familinx-raw/relations-anon.txt")) |> 
  group_by(child) |> 
  summarise(parents = str_c(parent, collapse = ",")) |>
  separate(parents, into = c("parent1", "parent2"), sep = ",") |>
  rename(profileid = child) |>
  as.data.table()

profile_subsets <- list.files(here("data/intermediate/"), full.names = TRUE) |> str_subset("[0-9]_assigned_new.csv")

profiles <- future_map_dfr(
  profile_subsets,
  function(x) {read_csv(x)}
)

df <- left_join(profiles, links) |>
  mutate(across(ends_with("place"), ~na_if(.x, ", ,")))

write_csv(df, file = here("data/intermediate/profiles_recombined_new.csv"))
