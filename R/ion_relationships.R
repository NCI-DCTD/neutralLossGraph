#' Calculate Pairwise Ion Relationships
#'
#' Examines m/z differences among features to infer potential charge states,
#' neutral losses, and adduct relationships.
#'
#' @param mzs Numeric vector of m/z values
#' @param rt Retention time (shared across all mzs in this set)
#' @param ids Vector of feature IDs corresponding to mzs
#' @param adduct_data Output from prepareAdductData(), table of common positive 
#'   adduct molecular weights and neutral losses for ms1 feature relationships
#'
#' @return Data.frame of inferred pairwise relationships
calculatePairWiseIons <- function(mzs, rt, ids, adduct_data) {
  adducts <- adduct_data$adducts
  adducts.pos <- adduct_data$adducts_pos
  adducts.tbl <- adduct_data$adducts_tbl

  pwIonTbl <- data.frame(matrix(ncol = 11, nrow = 0), stringsAsFactors = FALSE)
  colnames(pwIonTbl) <- c('rt', 'mz_1', 'mz_2', 'mz_diff', 'z', 'z_adduct',
                          'loss', 'ion_1', 'ion_2', 'name_1', 'name_2')

  for (m in 1:(length(mzs) - 1)) {
    for (n in (m + 1):length(mzs)) {
      zCalcs <- round(1 / sapply(adducts.pos$WT, function(x) (mzs[m] - x) / (mzs[n] - x)), 3)
      zIdx <- which(zCalcs %in% c(2, 3, 4))
      z <- zCalcs[zIdx][1]
      z_adduct <- adducts.pos$Adduct[zIdx][1]

      lossCalcs <- round(sapply(adducts[adducts$charge == 0, ]$WT, function(x) (mzs[n] - mzs[m]) / x), 3)
      lossIdx <- which(lossCalcs %in% c(1, 2, 3))
      loss <- if (length(lossIdx) > 0) paste0(lossCalcs[lossIdx][1], '*', adducts$Adduct[lossIdx][1]) else NA

      diffCalcs <- round(sapply(adducts.tbl$diff, function(x) (mzs[n] - mzs[m]) / x), 3)
      diffIdx <- which(diffCalcs == 1)
      ion_1 <- if (length(diffIdx)) adducts.tbl$adduct_2[diffIdx][1] else NA
      ion_2 <- if (length(diffIdx)) adducts.tbl$adduct_1[diffIdx][1] else NA

      pwIonTbl[nrow(pwIonTbl) + 1, ] <- c(rt, mzs[m], mzs[n], mzs[n] - mzs[m], z, z_adduct, loss, ion_1, ion_2, ids[m], ids[n])
    }
  }

  pwIonTbl
}

#' Annotate Features Based on Ion Relationships
#'
#' Adds adduct, neutral loss, and ion-pair annotations to features.
#'
#' @param tbl Feature table (with mz, id, etc.)
#' @param pw Pairwise ion table from calculatePairWiseIons
#' @return tbl with updated 'annotation' and 'direction'
annotateIonTbl <- function(tbl, pw) {

  for (i in seq_len(nrow(tbl))) {
    annotations <- pw[pw$mz_1 == tbl[i, 'mz'], ]

    for (a in seq_len(nrow(annotations))) {
      row <- annotations[a, ]
      idx_target <- which(tbl$mz == row$mz_2)

      if (!is.na(row$z)) {
        z <- as.numeric(row$z)
        charge <- tbl$charge[i]
        charge <- ifelse(!is.na(charge) && charge > 1, as.character(charge), "")
        # tbl$annotation[i] <- paste0(tbl$annotation[i], "[M + ", charge, row$z_adduct, "] ")
        tbl$annotation[i] <- paste0(tbl$annotation[i], row$name_1, ":[M + ", charge, row$z_adduct, "]->",row$name_2, ":[", z, "M + ", row$z_adduct, "],")
        tbl$direction[i] <- paste0(tbl$direction[i], row$name_1, "->", row$name_2, ",")
        if (length(idx_target) > 0) {
          tbl$annotation[idx_target] <- paste0(tbl$annotation[idx_target], "[", z, "M + ", row$z_adduct, "] ")
          # tbl$direction[idx_target] <- paste0(tbl$direction[idx_target], row$name_1, "<-", row$name_2, ",") # added for clarity
        }
      }

      if (!is.na(row$loss)) {
        tbl$annotation[i] <- paste0(tbl$annotation[i], "[", row$mz_2, " - ", row$loss, "] ")
        tbl$direction[i] <- paste0(tbl$direction[i], row$name_1, "->", row$name_2, ",")
      }

      if (!is.na(row$ion_1) && !is.na(row$ion_2)) {
        # tbl$annotation[i] <- paste0(tbl$annotation[i], "[M + ", row$ion_1, "] ")
        tbl$annotation[i] <- paste0(tbl$annotation[i], row$name_1, ":[M + ", row$ion_1, "]->",row$name_2, ":[M + ", row$ion_2, "],")
        if (length(idx_target) > 0) {
          tbl$annotation[idx_target] <- paste0(tbl$annotation[idx_target], "[M + ", row$ion_2, "] ")
          # tbl$direction[idx_target] <- paste0(tbl$direction[idx_target], row$name_1, "<-", row$name_2, ",") # added for clarity
        }
        tbl$direction[i] <- paste0(tbl$direction[i], row$name_1, "->", row$name_2, ",")
      }
    }
  }
  tbl[c("annotation", "direction")]
}
