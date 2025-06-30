#' Parse MGF-style text into a list of metadata and spectra blocks
#'
#' This function extracts and parses individual blocks of mass spectrometry data
#' enclosed between `BEGIN IONS` and `END IONS` tags in MGF-style text. Each block
#' is split into named metadata and a numeric matrix of m/z and intensity values.
#'
#' @param text A single character string containing the entire MGF-style content.
#'
#' @return A list of length-N, where each element is a list with:
#'   - `metadata`: a named character vector of metadata fields (e.g., FEATURE_ID, PEPMASS)
#'   - `spectrum`: a numeric matrix with two columns: mz and intensity
#'
#' @examples
#' mgf_text <- "BEGIN IONS\nFEATURE_ID=1\nPEPMASS=123.4\n10 100\n20 200\nEND IONS"
#' parse_mgf_blocks(mgf_text)
#'
#' @export
parse_mgf_blocks <- function(text) {
  # Split the text by BEGIN/END tags
  blocks <- unlist(strsplit(text, "BEGIN IONS\\s*|\\s*END IONS"))
  
  # Remove whitespace-only and empty strings
  blocks <- blocks[nzchar(trimws(blocks))]
  
  # Parse each block
  parsed <- lapply(blocks, function(block) {
    lines <- strsplit(block, "\n")[[1]]
    lines <- trimws(lines[nzchar(trimws(lines))])
    
    # Separate metadata and spectrum lines
    metadata_lines <- lines[grepl("=", lines)]
    spectrum_lines <- lines[!grepl("=", lines)]
    
    # Parse metadata into named character vector
    metadata <- setNames(
      lapply(strsplit(metadata_lines, "="), function(x) x[2]),
      sapply(strsplit(metadata_lines, "="), function(x) x[1])
    )
    
    # Parse spectrum lines into numeric matrix
    spectrum <- do.call(rbind, lapply(spectrum_lines, function(line) {
      as.numeric(strsplit(line, "\\s+")[[1]])
    }))
    
    if (!is.null(spectrum)) {
      colnames(spectrum) <- c("mz", "intensity")
    }
    
    list(metadata = metadata, spectrum = spectrum)
  })
  
  parsed
}
