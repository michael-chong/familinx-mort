library(tidyverse)
library(countrycode)

country_codes <- list()

## United States 
country_codes$USA <- countrycode::codelist |>
  filter(country.name.en == "United States") |>
  select(where(is.character)) |>
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |>
  filter(str_detect(name, "name")) |> 
  pull(value) |>
  str_to_upper() |> c(
  ", US$", 
  "USA", 
  "UNITED STATES", 
  "UNITED STATES OF AMERICA",  
  ", AMERICA", 
  'REPUBLIC OF AMERICA', 
  "NEW ENGLAND",
  "UNITED COLONIES OF AMERICA",
  "MASSACHUSETTS",
  "CONNECTICUT",
  "ILLINOIS",
  "NEW YORK",
  "NORTH CAROLINA",
  "SOUTH CAROLINA",
  "NEW HAMPSHIRE",
  "KENTUCKY",
  "RHODE ISLAND",
  "COLONIAL AMERICA",
  "PROVINCE OF VIRGINIA",
  "TEXAS", 
  "CALIFORNIA"
) |>
  unique()

## France
country_codes$FRA <- countrycode::codelist |>
  filter(country.name.en == "France") |>
  select(where(is.character)) |>
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |>
  filter(str_detect(name, "name")) |> 
  pull(value) |>
  str_to_upper() |> 
  str_subset("FRANCE", negate = TRUE) |>
  c(", FRA$", "FR$", "\\w*(?<!NEW )FRANCE", "\\w*(?<!NOUVELLE[ -])FRANCE", "R(É|E)P.*FRANÇAISE", "FRENCH.?REPUBLIC")

## Finland
country_codes$FIN <- countrycode::codelist |>
  filter(country.name.en == "Finland") |>
  select(where(is.character)) |>
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |>
  filter(str_detect(name, "name")) |> 
  pull(value) |>
  str_to_upper() |>
  c(" FIN$", "FI$","SOOMLANE")

## Denmark
country_codes$DNK <- countrycode::codelist |>
  filter(country.name.en == "Denmark") |>
  select(where(is.character)) |>
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |>
  filter(str_detect(name, "name")) |> 
  pull(value) |>
  str_to_upper() |>
  c("DK$", "DNK")

## Norway
country_codes$NOR <- countrycode::codelist |>
  filter(country.name.en == "Norway") |>
  select(where(is.character)) |>
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |>
  filter(str_detect(name, "name")) |> 
  pull(value) |>
  str_to_upper() |>
  c(" NOR$", " NO$")

## Sweden
country_codes$SWE <-  countrycode::codelist |>
  filter(country.name.en == "Sweden") |>
  select(where(is.character)) |>
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |>
  filter(str_detect(name, "name")) |> 
  pull(value) |>
  str_to_upper() |>
  c(" SWE$", " SE$", "SUÃ¨DE")

## Belgium
country_codes$BEL <- countrycode::codelist |> 
  filter(country.name.en == "Belgium") |> 
  select(where(is.character)) |> 
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |> 
  filter(str_detect(name, "name")) |>
  pull(value) |>
  str_to_upper() |>
  c(" BE$", " BEL$")

## Netherlands
country_codes$NLD <- countrycode::codelist |> 
  filter(country.name.en == "Netherlands") |> 
  select(where(is.character)) |> 
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |> 
  filter(str_detect(name, "name")) |>
  pull(value) |>
  str_to_upper() |>
  str_subset("^OLAND$|^HOLLAND*", negate = TRUE) |>  # remove these because 
  c(" NL$", " NLD$", "\\w*(?<!NEW )HOLLAND")

## Great Britain
country_codes$GBR <- countrycode::codelist |> 
  filter(country.name.en == "United Kingdom") |> 
  select(where(is.character)) |> 
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |> 
  filter(str_detect(name, "name")) |>
  pull(value) |>
  str_to_upper() |>
  str_subset("^ENGLAND$", negate = TRUE) |> 
  c("BRITAIN", "\\w*(?<!NEW )SCOTLAND", " UK", " GBR", "\\w*(?<!NEW )ENGLAND", "\\w*(?<!NEW SOUTH )WALES")

## Switzerland
country_codes$CHE <- countrycode::codelist |>
  filter(country.name.en == "Switzerland") |>
  select(where(is.character)) |> 
  pivot_longer(!country.name.en) |>
  distinct(value, .keep_all = TRUE) |> 
  filter(str_detect(name, "name")) |>
  pull(value) |>
  str_to_upper() |>
  c(" CH$", " CHE$")

write_rds(country_codes, here::here("data/country_codes.rds"))
