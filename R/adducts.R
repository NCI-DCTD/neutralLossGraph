#' Prepare Adduct Data for Annotation
#'
#' Loads adduct information and calculates pairwise mass differences
#' for positive-mode ions to support downstream annotation.
#'
#' @param adducts_csv_path Path to a CSV file with columns: Adduct, WT, charge
#' @return A list containing:
#'   - `adducts`: all adducts
#'   - `adducts_pos`: positive-mode adducts
#'   - `adducts_tbl`: pairwise mass differences between positive adducts
prepareAdductData <- function(adducts_csv_path = system.file("extdata", "adducts.csv", package = "msannotator")) {
  adducts <- read.csv(adducts_csv_path, stringsAsFactors = FALSE)
  adducts.pos <- adducts[adducts$charge == 1, ]
  
  adducts.tbl <- data.frame(adduct_1 = character(), adduct_2 = character(), diff = numeric(), stringsAsFactors = FALSE)
  
  for (m in seq(nrow(adducts.pos), 2)) {
    for (n in seq(m - 1, 1)) {
      adducts.tbl[nrow(adducts.tbl) + 1, ] <- list(
        adducts.pos[m, 'Adduct'],
        adducts.pos[n, 'Adduct'],
        adducts.pos[m, 'WT'] - adducts.pos[n, 'WT']
      )
    }
  }
  
  list(
    adducts = adducts,
    adducts_pos = adducts.pos,
    adducts_tbl = adducts.tbl
  )
}
