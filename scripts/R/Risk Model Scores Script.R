# ============================================================
# ETHE1 Revised Variant Prioritization Risk Model
# ============================================================
#
# Input:
#   results/prioritization/ETHE1 Variants Risk Model Results.csv
#
# Output:
#   results/prioritization/ETHE1_NEW_MODEL.csv
#
# Dependencies:
#   readr, dplyr
#
# ============================================================

library(dplyr)
library(readr)

# ------------------------------------------------------------
# 1. Define input and output paths
# ------------------------------------------------------------

input_file <- file.path(
  "results",
  "prioritization",
  "ETHE1 Variants Risk Model Results.csv"
)

output_file <- file.path(
  "results",
  "prioritization",
  "ETHE1_NEW_MODEL.csv"
)

# ------------------------------------------------------------
# 2. Validate input file
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

ETHE1 <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "Consequence",
  "Active_Site",
  "Iron_Binding",
  "Epitope",
  "CADD_PHRED",
  "Position"
)

missing_columns <- setdiff(
  required_columns,
  names(ETHE1)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 3. Calculate revised prioritization scores
# ------------------------------------------------------------

ETHE1_Revised <- ETHE1 %>%
  mutate(

    # Consequence score
    Consequence_Score = case_when(
      Consequence == "frameshift_variant" ~ 5,
      Consequence == "stop_gained" ~ 5,
      Consequence == "start_lost" ~ 4,
      Consequence == "stop_lost" ~ 4,
      Consequence %in% c(
        "inframe_insertion",
        "inframe_deletion"
      ) ~ 4,
      Consequence == "missense_variant" ~ 2,
      Consequence == "synonymous_variant" ~ 0,
      TRUE ~ 0
    ),

    # Functional annotation scores
    Active_Site_Score = if_else(
      coalesce(Active_Site, FALSE),
      4,
      0
    ),

    Iron_Binding_Score = if_else(
      coalesce(Iron_Binding, FALSE),
      4,
      0
    ),

    Epitope_Score = if_else(
      coalesce(Epitope, FALSE),
      1,
      0
    ),

    # CADD score
    CADD_Score = case_when(
      CADD_PHRED >= 29 ~ 5,
      CADD_PHRED >= 20 ~ 4,
      CADD_PHRED >= 10 ~ 3,
      CADD_PHRED < 10 ~ 0,
      TRUE ~ 0
    ),

    # Initiator methionine score
    M1_Score = if_else(
      Position == 1,
      3,
      0
    ),

    # Second-order structural proximity score
    Second_Order_Score = case_when(

      Position %in% c(
        33, 34, 35, 78, 81, 84, 136, 139,
        165, 166, 174, 193, 194, 195, 214,
        215, 222, 226, 231, 234
      ) ~ 2,

      Position %in% c(
        52, 80, 82, 106, 152, 155, 156, 158,
        164, 167, 196, 198, 212, 229, 232,
        233, 236
      ) ~ 1,

      TRUE ~ 0
    )
  ) %>%

  # ----------------------------------------------------------
  # 4. Calculate total revised priority score
  # ----------------------------------------------------------

  mutate(
    Priority_Score_Revised =
      Consequence_Score +
      Active_Site_Score +
      Iron_Binding_Score +
      Epitope_Score +
      CADD_Score +
      M1_Score +
      Second_Order_Score
  ) %>%

  # ----------------------------------------------------------
  # 5. Assign revised priority categories
  # ----------------------------------------------------------

  mutate(
    Priority_Revised = case_when(
      Priority_Score_Revised >= 12 ~ "Very High",
      Priority_Score_Revised >= 8 ~ "High",
      Priority_Score_Revised >= 5 ~ "Moderate",
      TRUE ~ "Low"
    )
  )

# ------------------------------------------------------------
# 6. Export results
# ------------------------------------------------------------

write_csv(
  ETHE1_Revised,
  output_file
)

# ------------------------------------------------------------
# 7. Print run summary
# ------------------------------------------------------------

cat("\nETHE1 revised risk model completed successfully.\n")
cat("Input variants:", nrow(ETHE1_Revised), "\n")
cat("Output file:", output_file, "\n\n")

cat("Priority distribution:\n")

print(
  ETHE1_Revised %>%
    count(Priority_Revised, name = "Variant_Count")
)

cat("\nConsequence score distribution:\n")

print(
  ETHE1_Revised %>%
    count(Consequence, Consequence_Score)
)
