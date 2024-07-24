library(tidyverse)
library(here)
library(mgcv)
library(cmdstanr)

source(here("code/functions.R"))

year_min <- 1800
year_max <- 1900
year_interval_length <- 5
year_breaks <- seq(year_min, year_max, by = year_interval_length)
n_yb <- length(year_breaks)

age_min <- 5
age_max <- 90
age_interval_length <- 5
age_breaks <- c(0, 1, seq(age_min, age_max, age_interval_length), Inf)
n_ab <- length(age_breaks)

country_ref <- c("FRA", "SWE", "DNK", "NOR", "BEL", "GBR", "NLD", "FIN", "CHE")
country_unk <- c("USA", "FRA")

deaths_hmd <- read_csv(here("data/clean/deaths_hmd_all.csv")) |>
  filter(country %in% country_ref) |>
  make_year_age_cells(
    year_breaks = year_breaks, 
    age_breaks = age_breaks, 
    out = "df", 
    gender_group = TRUE
  ) |>
  filter(!is.na(year)) |> 
  rename(deaths = n)

exposure_hmd <- read_csv(here("data/clean/exposure_hmd_all.csv")) |>
  filter(country %in% country_ref) |>
  make_year_age_cells(
    year_breaks = year_breaks, 
    age_breaks = age_breaks, 
    out = "df", 
    gender_group = TRUE
  ) |>
  filter(!is.na(year)) |> 
  rename(exposure = n)

df_hmd <- left_join(deaths_hmd, exposure_hmd, relationship = "one-to-one") |>
  filter(!(country == "FRA"))


deaths_fam <- read_rds(here("data/clean/migrant_deaths_grouped_new.rds")) |>
  filter(!is.na(age), !is.na(country), !is.na(year)) |> 
  filter(as.numeric(as.character(year)) >= year_min) |> 
  filter(country %in% c(country_ref, country_unk)) |>
  group_by(country, age_group = age, year, gender) |>
  summarise(deaths = sum(n))

exposure_fam <- read_rds(here("data/clean/migrant_exposures_grouped_new.rds")) |>
  filter(!is.na(age), !is.na(country), !is.na(year)) |>
  filter(as.numeric(as.character(year)) >= year_min) |>
  filter(country %in% c(country_ref, country_unk)) |>
  group_by(country, age_group = age, year, gender) |>
  summarise(exposure = sum(n))

df_fam <- left_join(deaths_fam, exposure_fam) |>
  arrange(country, gender, year, age_group) |>
  mutate(country = as.factor(country)) |> 
  mutate(row_id = row_number()) |> 
  group_by(country, gender) |>
  mutate(group_id = cur_group_id()) |>
  ungroup()

group_gender <- df_fam |>
  arrange(group_id) |> 
  distinct(group_id, gender) |>
  pull(gender)

country_group <- df_fam |>
  arrange(group_id) |>
  distinct(group_id, country) |>
  pull(country)

# replicates roughly what is done under the hood in brms GAMs
# generate spline terms that can be estimated as random effects
## decomposes into fixed linear trends and random effect components
# see mgcv::t2() to choose a different kind of basis/penalty
## this uses the default (cubic regression splines, penalized by the "integrated square 2nd derivative")
sm_psi <- smoothCon(
  t2(as.numeric(age_group), as.numeric(year), k = c(10, 9)),
  df_fam, 
  knots = list(
    "as.numeric(age_group)" = c(2, 2.5, 3, 3.5, 4, 6, 8, 12, 16, 18)
  ),
  absorb.cons = TRUE, 
  modCon = 3, 
  diagonal.penalty = TRUE
)[[1]]  |> mgcv::smooth2random(vnames = "", type = 2)

X <- sm_psi$Xf
Z <- sm_psi$rand

sm_pc <- smoothCon(
  s(as.numeric(year), k = 4),
  df_fam |> distinct(country, gender, year) |> arrange(country, gender, year), 
  absorb.cons = TRUE, 
  modCon = 3, 
  diagonal.penalty = TRUE
)[[1]]  |> mgcv::smooth2random(vnames = "", type = 2)

X_pc <- sm_pc$Xf
Z_pc <- sm_pc$rand

mat_fam_deaths <- df_fam |> reshape2::acast(country*gender*year ~ age_group, value.var = "deaths")
mat_fam_exposure <- df_fam |> reshape2::acast(country*gender*year ~ age_group, value.var = "exposure")

year_fam <- match(
  str_extract(rownames(mat_fam_deaths), str_flatten(year_breaks, "|")),
  year_breaks
)

df_hmd <- df_hmd |>
  left_join(
    df_fam |> select(-deaths, -exposure),
    by = c("country", "gender", "year", "age_group"),
    relationship = "one-to-one"
  ) 

mat_hmd_deaths <- df_hmd |> reshape2::acast(country*gender*year ~ age_group, value.var = "deaths")
mat_hmd_exposure <- df_hmd |> reshape2::acast(country*gender*year ~ age_group, value.var = "exposure")

ii <- match(rownames(mat_hmd_deaths), rownames(mat_fam_deaths))

# country-sex groups in long format
groups_mat <- df_fam |> reshape2::acast(country*gender*year ~ age_group, value.var = "group_id")

hmd_pc_male <- log(mat_hmd_deaths / mat_hmd_exposure) |>
  (\(.) .[str_detect(rownames(mat_hmd_deaths), "_male_"), ])() |> 
  svd() |>
  with(v) |>
  (\(.) .[,1:4] )() |>
  t()

hmd_pc_female <- log(mat_hmd_deaths / mat_hmd_exposure) |>
  (\(.) .[str_detect(rownames(mat_hmd_deaths), "_female_"), ])() |> 
  svd() |>
  with(v) |>
  (\(.) .[,1:4] )() |>
  t()

# extra objects to produce life expectancy estimates for the unknown country
j <- which(str_detect(rownames(mat_fam_deaths), str_flatten(country_unk, "|")))

invalid <- (rownames(mat_fam_deaths) |> str_detect("NLD")) & (year_fam <= 4)
valid_rows <- which(!invalid)

data_list <- list(
  n_countries = n_distinct(c(country_ref, country_unk)),
  n_periods = length(year_breaks) - 1,
  n_ages = length(age_breaks) - 1,
  N_fam = nrow(mat_fam_deaths),
  N_long = nrow(mat_fam_deaths)*(length(age_breaks)-1),
  d_fam = mat_fam_deaths,
  P_fam = mat_fam_exposure,
  n_valid = length(valid_rows),
  valid_rows = valid_rows,
  year_fam = year_fam,
  country_fam = str_extract(rownames(mat_fam_deaths), str_c(c(country_ref, country_unk), collapse = "|")) |> match(c(country_ref, country_unk)),
  country_group = country_group |> match(c(country_ref, country_unk)),
  gender_fam = 1*str_detect(rownames(mat_fam_deaths), "_male_") + 1,
  gender_long = 1*(df_fam$gender == "male") + 1,
  gender_group = 1*(group_gender == "male") + 1,
  male_long_ind  = which(df_fam$gender == "male"),
  female_long_ind = which(df_fam$gender == "female"),
  n_groups = max(df_fam$group_id), 
  group_long = df_fam$group_id,
  group_ind = sapply(1:max(df_fam$group_id), function(x) which(df_fam$group_id == x)) |> t(),
  group_fam = groups_mat[,1],
  N_hmd = nrow(mat_hmd_deaths),
  d_hmd = mat_hmd_deaths,
  P_hmd = mat_hmd_exposure,
  country_hmd = match(df_hmd$country, c(country_ref, country_unk)),
  ii = ii,
  # for life expectancy calculations
  n_agegap = lead(age_breaks[-n_ab], default = 999999) - age_breaks[-n_ab],
  n_groups_unk = length(j),
  j = j,
  # spline stuff for adjustment factor
  k_X =  dim(X)[2], 
  k_Z1 = dim(Z[[1]])[2],
  k_Z2 = dim(Z[[2]])[2],
  k_Z3 = dim(Z[[3]])[2],
  X = X,
  Z1 = Z[[1]],
  Z2 = Z[[2]],
  Z3 = Z[[3]],
  scale_random_spline = 1,
  # svd stuff for mortality rate
  k_pc = nrow(hmd_pc_male),
  V_pc_female = hmd_pc_female,
  V_pc_male = hmd_pc_male,
  k_X_pc = ncol(X_pc),
  X_pc = X_pc,
  k_Z_pc = ncol(Z_pc[[1]]),
  Z_pc = Z_pc[[1]],
  # horseshoe
  hs_df = 1,
  hs_df_global = 1,
  hs_df_slab = 1,
  hs_scale_global = 0.01, 
  hs_scale_slab = 1,
  # other
  phi_numer_hmd = 1000,
  phi_numer_fam = 100
)

attr(data_list, "indices")  <- list(
  countries = c(country_ref, country_unk),
  periods = year_breaks,
  ages = age_breaks
)

attr(data_list, "dataframes") <- list(
  fam = df_fam,
  hmd = df_hmd
)

# recompile because potentially running on different computers/clusters; may accidentally transfer executable
mod <- cmdstan_model(here("code/modelling/model.stan"), pedantic = TRUE)

print(strrep("----", 10))
print(mod$model_name())
print(strrep("----", 10))

fit <- mod$sample(
  data = data_list,
  seed = 96,
  parallel_chains = 4,
  max_treedepth = 10,
  output_dir = here("output/draws/"),
  output_basename = "model_fit_i"
)

summary_table <- fit$summary()

write_rds(
  list(
    fit = fit,
    summary_table = summary_table,
    data = data_list
  ),
  file = here("output/intermediate/model_fit_i.rds")
)
