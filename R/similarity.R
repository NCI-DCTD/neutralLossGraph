#' Calculate Pairwise Fragment Similarity Between Two Spectra
#'
#' This function compares two MS/MS spectra by evaluating the similarity
#' between their internal fragment m/z difference patterns (i.e., pairwise distances).
#' It returns a score between 0 and 1 representing how similar the fragmentation
#' patterns are based on their internal structure.
#'
#' The method computes pairwise distances within each spectrum, filters out distances
#' less than or equal to 10 (to reduce noise), and compares the sets of distances.
#' The output is the higher of two metrics:
#' \enumerate{
#'   \item The proportion of total distance overlap.
#'   \item The proportion of distance values that match.
#' }
#'
#' @param v1 Numeric vector. m/z values of one MS/MS spectrum (e.g., fragment peaks).
#' @param v2 Numeric vector. m/z values of a second MS/MS spectrum to compare against.
#'
#' @return A numeric similarity score between 0 (no similarity) and 1 (perfect match).
#'
#' @details
#' The shorter of the two input vectors is used as the reference. Pairwise distances
#' are computed using `dist()`, rounded to 2 decimal places, and filtered to retain
#' only those >10.
#'
#' If either vector contains the exact value `10`, pairwise distance calculation
#' is skipped, and that vector is used as-is. Useful for adding known neutral loss vectors
#'
#' @examples
#' calculatePWDiff(c(100, 150, 200), c(100, 150, 250))
#' calculatePWDiff(c(101, 111, 121, 131), c(101, 112, 123, 134))
#'
#' @export
calculatePWDiff <- function(v1, v2) {
  len1 <- length(v1)
  len2 <- length(v2)

  # Assign shorter vector to A, longer to C
  if (len1 <= len2) {
    A <- v1
    C <- v2
  } else {
    A <- v2
    C <- v1
  }

  # Compute pairwise distances or handle edge case with '10'
  A_dists <- if (10 %in% A) A else dist(A)
  C_dists <- if (10 %in% C) C else dist(C)

  # Round and filter
  A_dists <- unique(round(A_dists[A_dists > 10], 2))
  C_dists <- unique(round(C_dists[C_dists > 10], 2))

  # Compute similarity score
  if (length(A_dists) > 0) {
    overlap_sum <- sum(A_dists[A_dists %in% C_dists])
    similarity <- max(
      overlap_sum / sum(A_dists),
      sum(A_dists %in% C_dists) / length(A_dists)
    )
    return(similarity)
  } else {
    return(0)
  }
}
