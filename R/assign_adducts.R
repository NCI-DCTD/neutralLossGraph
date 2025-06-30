#' Assign Adduct Annotations to Mass Spectrometry Features
#'
#' Clusters MS features by retention time and m/z, identifies related ions, and
#' annotates them with likely adducts or relationships. Also estimates neutral masses
#' based on annotation patterns and configurable correction rules.
#'
#' @param rt A data.frame with MS feature information. Must include:
#'   - `mz`: mass-to-charge ratio
#'   - `rt`: retention time
#'   - `id`: a unique identifier
#' @param rt_h Numeric. Clustering height for retention time grouping. Default = 30.
#' @param mz_h Numeric. Clustering height for m/z grouping. Default = 0.005.
#' @param neutral_rules Named numeric vector. Maps adduct labels to neutral mass
#'   correction values. Default includes H, Na, NH4, and 2M+H (dimer).
#' @param includeMS1 Logical. If TRUE, directional MS1-based annotations are calculated
#'   based on retention time + mass to charge group and known adduct values. Default is TRUE.
#'
#' @return A data.frame with added columns:
#'   - `rtGroup`, `mzGroup`: cluster assignments
#'   - `annotation`, `direction`: inferred adduct info
#'   - `Neutral`: estimated neutral mass
#'
#' @export

assignAdducts <- function(rt, rt_h = 30, mz_h = 0.005,
                          neutral_rules = c(
                            "M + H"   = 1.007276,
                            "M + Na"  = 22.989218,
                            "M + NH4" = 18.033823,
                            "2M + H"  = 1.007276  # used when charge is NA
                          ), includeMS1 = TRUE) {

  rt$rtGroup <- 1
  rt$mzGroup <- 1
  rt$origOrder <- seq_len(nrow(rt))


  if (nrow(rt) > 1) {
    # RT clustering
    rt <- rt[order(rt$rt), ]
    rt_cuts <- cutree(hclust(dist(rt$rt), method = "centroid"), h = rt_h)
    if (length(rt_cuts) == nrow(rt)) rt$rtGroup <- rt_cuts

    # m/z clustering
    rt <- rt[order(rt$mz), ]
    mz_cuts <- cutree(hclust(dist(rt$mz), method = "centroid"), h = mz_h)
    if (length(mz_cuts) == nrow(rt)) rt$mzGroup <- mz_cuts
  }

  rt <- rt[order(rt$origOrder), ]
  rt$annotation <- ""
  rt$direction <- ""

  if (includeMS1) {
    adduct_data <- prepareAdductData()
    # -------------------------
    # Annotate adducts
    # -------------------------
    for (i in unique(rt$rtGroup)) {
      rtSet <- which(rt$rtGroup == i)
      if (length(rtSet) > 1) {
        pw <- calculatePairWiseIons(rt[rtSet, ]$mz, i, rt[rtSet, ]$id, adduct_data)

        annotations <- annotateIonTbl(rt[rtSet, ], pw)
        rt[rtSet, "annotation"] <- stringr::str_trim(annotations$annotation)
        rt[rtSet, "direction"] <- stringr::str_trim(annotations$direction)
      } else {
        rt[rtSet, "annotation"] <- "[M + H]"
      }
    }
    # -------------------------
    # Estimate Neutral Mass
    # -------------------------
    rt$Neutral <- NA
    for (pattern in names(neutral_rules)) {
      matches <- grepl(paste0("[", pattern, "]"), rt$annotation, fixed = TRUE)

      # Handle special case for '2M + H' (only when charge is NA)
      if (pattern == "2M + H") {
        matches <- matches & is.na(rt$charge)
      }

      rt[matches, "Neutral"] <- rt[matches, "mz"] - neutral_rules[[pattern]]
    }
  }
  return(rt)
}
