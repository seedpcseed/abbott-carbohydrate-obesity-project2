suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(stringr)
})

micro_path <- "../microbiome_analysis/zr24558.16S_250813.zymo/00...AllSamples.Bac16Sv34/Sample_Information/sample.metadata.csv"
cat("reading:", micro_path, "\n")
micro_metadata_raw <- readr::read_tsv(micro_path, show_col_types = FALSE) %>%
  janitor::clean_names()

micro_metadata <- micro_metadata_raw %>%
  dplyr::mutate(
    customer_label = stringr::str_replace_all(customer_label, "\\.\\.", "."),
    customer_label = stringr::str_trim(customer_label)
  ) %>%
  tidyr::separate(
    customer_label,
    into = c("aliquot_id", "replicate_code", "timepoint"),
    sep = "\\.",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    aliquot_id = stringr::str_to_upper(aliquot_id),
    donor_id = stringr::str_extract(aliquot_id, "^[0-9]+"),
    subject_id = aliquot_id
  )

odd <- micro_metadata %>%
  dplyr::distinct(customer_label, aliquot_id, donor_id) %>%
  dplyr::filter(
    is.na(donor_id) |
      !stringr::str_detect(aliquot_id, "^[0-9]+[A-Z]$")
  )
print(odd)

study <- micro_metadata %>%
  dplyr::filter(stringr::str_detect(aliquot_id, "^[0-9]+[A-Z]$"))

cat(
  "study donors", dplyr::n_distinct(study$donor_id),
  "study aliquots", dplyr::n_distinct(study$aliquot_id), "\n"
)
stopifnot(dplyr::n_distinct(study$donor_id) == 16)
stopifnot(dplyr::n_distinct(study$aliquot_id) == 32)
row <- study %>%
  dplyr::filter(customer_label == "84B.R1.48H") %>%
  dplyr::slice(1)
stopifnot(as.character(row$donor_id) == "84", as.character(row$aliquot_id) == "84B")
cat(
  "smoke OK:", as.character(row$customer_label), "->",
  as.character(row$donor_id), as.character(row$aliquot_id), "\n"
)
