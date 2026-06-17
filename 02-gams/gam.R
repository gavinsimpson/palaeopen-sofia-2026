# scratch code script for preparing data and exploring

pkgs <- c(
  "here", "mgcv", "readr", "gratia", "tidyr"
)

vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)

## code for swap chemistry
swapc <- readCEP("~/work/data/swap/SWAP178.ENV", long = TRUE)

swapc |>
  as_tibble() |>
  mutate(
    row = stringr::str_replace(row, "^X", ""),
    col = replace_values(
      col,
      "Ph" ~ "pH"
    )
  ) |>
  rename(
    sample_id = row,
    variable = col,
    value = val
  ) |>
  pivot_wider(
    id_cols = "sample_id",
    names_from = "variable",
    values_from = "value"
  ) |>
  write_csv(file = here("00-data/swap-chemistry.csv"))

## code for swap diatoms
swapd <- read_csv2(here("00-data/swap-diatom-counts.csv")) |>
  janitor::clean_names() |>
  select(
    c(original_diatom_sample_id, taxon_id, abun)
  ) |>
  rename(
    sample_id = original_diatom_sample_id
  ) |>
  pivot_wider(
    id_cols = "sample_id",
    names_from = "taxon_id",
    values_from = "abun"
  ) |>
  mutate(
    across(
      AC013A:NA057A,
      .fns = ~ replace_na(.x, replace = 0)
    )
  ) |>
  rowwise() |>
  mutate(
    total_count = sum(c_across(AC013A:NA057A))
  )

# code for processed swap chemistry
swapc <- read_csv(here("00-data/swap-chemistry.csv"))

swap_sub <- swapd |>
  select(
    c("sample_id", "AC013A", "BR006A", "FR005D", "total_count")
  ) |>
  left_join(
    swapc |> select(c("sample_id", "pH", "Altot", "DOC", "Alkal")),
    by = join_by("sample_id")
  )

m_ph <- gam(
  cbind(FR005D, total_count - FR005D) ~ s(pH, k = 6),
  data = swap_sub,
  family = binomial(),
  method = "REML"
)

m_ph |> draw()

m_ph |> conditional_values(
  condition = "pH"
) |>
draw()

m_doc <- gam(
  cbind(FR005D, total_count - FR005D) ~ s(DOC, k = 6),
  data = swap_sub,
  family = binomial(),
  method = "REML"
)

m_doc |> draw()

m_doc |> conditional_values(
  condition = "DOC"
) |>
draw()
