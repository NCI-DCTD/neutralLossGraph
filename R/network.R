#' Run a test of neutral loss analysis using example MGF files
#'
#' This function locates `.mgf` files in the `extdata` folder of the `mypackage`
#' package and runs them through `NeutralLossGraph::networkMGFs()` to generate a network.
#'
#' @details
#' This is a helper function intended for testing or demonstrating the
#' `networkMGFs` function in the `NeutralLossGraph` package. The generated GML
#' output is written to `package_test` in current working directory.
#'
#' @return No return value. Called for side effect of writing a GML file.
#' @export
#'
#' @examples
#' graphMinalemine()
graphMinalemine <- function() {
  # Get the path to the extdata folder of the installed package
  extdata_path <- system.file("extdata", package = "NeutralLossGraph")

  # List all .mgf files
  mgf_files <- list.files(extdata_path, pattern = "\\.mgf$", full.names = TRUE)

  # Run network generation
  NeutralLossGraph::networkMGFs(
    mgf_files,
    outputGML = 'package_test',
    cutoff = 0.5,
    includeFragments = FALSE,
    includeMS1 = FALSE
  )
}

#' Generate a Molecular Similarity Network from MGF Files
#'
#' This function processes one or more MGF files and builds a molecular network based on
#' precursor and fragment mass spectra similarity. It filters, scores, annotates, and clusters
#' MS/MS spectra, then exports the result as a GraphML file for visualization.
#'
#' @param files A character vector of paths to MGF files.
#' @param cutoff Numeric threshold for pairwise similarity. Only spectra with a similarity
#'   above this threshold are included in the network. Default is 0.5.
#' @param outputGML The file path prefix for saving the output GraphML file. Default is
#'   `'msms_network'` in current working directory.
#' @param includeMS1 Logical. If TRUE, directional MS1-based annotations are used to modify
#'   edge weights in the network. Annotations will be added to network edges. Default is FALSE.
#' @param includeFragments Logical. If TRUE, pairwise calculations of two features will 
#'   include fragment values as well as neutral losses. Default is FALSE.
#' @param rt_h Numeric. Clustering height for retention time grouping. Default = 30.
#' @param mz_h Numeric. Clustering height for m/z grouping. Default = 0.005.
#'
#' @return No return value. Writes a GraphML file to the specified path.
#'
#' @details
#' This function computes pairwise similarity between MS/MS spectra based on fragment
#' mass-to-charge (m/z) values. It also incorporates MS1 annotation information if
#' specified, modifying edge weights for more meaningful clustering. The final result
#' is a similarity network stored in the GraphML format, which can be visualized using
#' tools such as Cytoscape.
#'
#' @examples
#' # Specify paths to MGF files and generate a similarity network
#' files <- list.files("data/mgf/", pattern = "\\.mgf$", full.names = TRUE)
#' networkMGFs(files, cutoff = 0.6, outputGML = "results/network", includeMS1 = TRUE)
#'
#' @import doParallel
#' @import foreach
#' @export
networkMGFs <- function(files, cutoff = 0.5, outputGML = 'msms_network', includeMS1 = FALSE, includeFragments = FALSE, rt_h = 30, mz_h = 0.005) {
  data     <- list()
  names    <- c()
  pmz      <- c()
  source   <- c()
  ioncount <- c()
  rt       <- c()
  charge   <- c()
  labels   <- c()

  # Register parallel backend
  n_cores <- parallel::detectCores() - 1
  cl      <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # Export necessary functions to workers
  parallel::clusterExport(cl, list("calculatePWDiff", "parse_mgf_blocks"))

  # Read and process MGF files in parallel
  results <- foreach(file = files, .packages = c("tools")) %dopar% {
    fs <- file.size(file)

    if (!is.na(fs) && fs > 0) {
      mgfdata <- parse_mgf_blocks(readChar(file, fs))
      scans   <- list()

      for (scan in mgfdata) {
        metadata <- scan$metadata
        pepmass  <- as.numeric(metadata$PEPMASS)
        spectra  <- scan$spectrum

        # remove product ions with intensity lower than 10 or are
        # greater mass than the parent ion
        filter  <- (spectra[, 'intensity'] > 10 & spectra[, 'mz'] <= pepmass)
        spectra <- spectra[filter,]

        if(is.null(dim(spectra))){
          mz <- spectra[1]
        } else {
            # remove any NA mz values. sort by decreasing intensity and select the top 30 remaining values         
            mz <- head(spectra[order(spectra[, 'intensity'], decreasing = TRUE), 'mz'], 30)
        }

        # final mz vector must contain at least 3 values
        if (length(mz) >= 3) {
          base_name <- tools::file_path_sans_ext(basename(file))

          # not all fields may be extracted from the mgf, but all that are
          # extracted are as.characters
          label     <- ifelse(is.null(metadata$NAME), NA, metadata$NAME)
          rt_val    <- ifelse(is.null(metadata$RTINSECONDS), 0, as.numeric(metadata$RTINSECONDS))
          scan_id   <- ifelse(is.null(metadata$FEATURE_ID), "0", metadata$FEATURE_ID)

          scans[[length(scans) + 1]] <- list(
            mz          = mz,
            precursorMz = pepmass,
            base_name   = base_name,
            scan_id     = scan_id,
            rt          = rt_val,
            charge      = metadata$CHARGE,
            label       = label
          )
        }
      }

      return(scans)
    } else {
      return(NULL)
    }
  }

  if (length(results) == 0) {
    warning("No data found in files. Exiting.")
    parallel::stopCluster(cl)
    return(NULL)
  }

  # Flatten the nested list of scans
  flat_results <- unlist(results, recursive = FALSE)

  for (res in flat_results) {
    data     <- append(data, list(c(res$mz, res$precursorMz)))
    pmz      <- c(pmz, res$precursorMz)
    names    <- c(names, paste(res$base_name, res$scan_id, sep = "_"))
    source   <- c(source, res$base_name)
    ioncount <- c(ioncount, length(res$mz))
    rt       <- c(rt, res$rt)
    charge   <- c(charge, res$charge)
    labels   <- c(labels, res$label)
  }

  # add zero to mz list for a scan to include fragment values for comparison
  if (includeFragments) {
    data <- lapply(data, function(x) c(x, 0))
  }

  if (length(pmz) == 0) {
    warning("Data of zero length. Exiting.")
    parallel::stopCluster(cl)
    return(NULL)
  }

  # Create similarity matrix
  n        <- length(data)
  pair_idx <- combn(n, 2, simplify = FALSE)

  sim_results <- foreach(pair = pair_idx, .combine = rbind) %dopar% {
    i   <- pair[1]
    j   <- pair[2]
    # a complete match (diff = 1) will be zeroed out/removed when 
    # constructing edge weights. set to a nominally lower
    # value to avoid losing the edge value
    sim <- min(calculatePWDiff(data[[i]], data[[j]]), 0.999)
    c(i = i, j = j, sim = sim)
  }

  # Stop cluster after similarity computation
  parallel::stopCluster(cl)

  tbl <- matrix(0, n, n, dimnames = list(names, names))

  for (k in 1:nrow(sim_results)) {
    i   <- sim_results[k, "i"]
    j   <- sim_results[k, "j"]
    sim <- sim_results[k, "sim"]
    tbl[i, j] <- sim
    tbl[j, i] <- sim
  }

  # Assign adduct annotations
  annotations <- assignAdducts(data.frame(
    mz     = pmz,
    rt     = rt,
    charge = charge,
    id     = names
    ), 
    includeMS1 = includeMS1,
    rt_h       = rt_h, 
    mz_h       = mz_h)

  # Optional MS1 edge boosting.
  # Extract all non-empty directions, split by comma
  aEdges <- unlist(strsplit(annotations$direction[annotations$direction != ""], ","), use.names = FALSE)

  # Loop over each direction and update the matrix. ignored if includeMS1 = FALSE
  for (edge in aEdges) {
    nodes <- strsplit(edge, "->", fixed = TRUE)[[1]]
    
    if (length(nodes) == 2) {
      from <- nodes[1]
      to   <- nodes[2]
      
      # Defensive check if names exist in the table
      if (from %in% rownames(tbl) && to %in% colnames(tbl)) {
        # set to the cutoff value to not supersede ms2 weights
        tbl[from, to] <- cutoff
      } else {
        warning(sprintf("Invalid edge: %s -> %s (not found in tbl)", from, to))
      }
    }
  }

  diag(tbl) <- 1
  # baseline cutoff for edges
  tbl[tbl < cutoff] <- 0

  # Create contraction index for vertex grouping
  rtMz <- unique(annotations[, c('rtGroup', 'mzGroup')])

  # Use match to get indices directly
  contract <- match(
    paste(annotations$rtGroup, annotations$mzGroup),
    paste(rtMz$rtGroup, rtMz$mzGroup)
  )

  # Write GraphML output
  outputGraphML(
    tbl,
    pmz,
    source,
    ioncount,
    annotations$annotation[annotations$annotation != ""],
    rt,
    contract,
    graph_name = outputGML,
    config = paste0('{ "cutoff" : ', cutoff, ', "includeMS1" : ', includeMS1, ', "includeFragments" : ', includeFragments, ' }' ),
    labels
  )
}


#' Export a Similarity Network to GraphML Format
#'
#' This function takes a similarity matrix, spectral annotations, and metadata to
#' generate a molecular similarity network as a GraphML file, suitable for visualization
#' in tools like Cytoscape. The network is based on pairwise fragment similarity between
#' MS/MS spectra.
#'
#' @param tbl A numeric matrix representing pairwise similarity scores between spectra.
#'   Similarity values should be between 0 (no similarity) and 1 (identical spectra).
#' @param pmzs A vector of precursor m/z values corresponding to the spectra.
#' @param source A vector containing the source file names for each spectrum.
#' @param ioncount A vector representing the number of ions (peaks) in each spectrum.
#' @param annotation A vector of annotation strings for each spectrum.
#' @param direction A vector describing the ion relationships for each spectrum.
#' @param rt A vector of retention times corresponding to each spectrum.
#' @param contract A vector indicating which vertices should be merged based on
#'   equivalence relations.
#' @param graph_name The desired name for the output GraphML file. Default is `'new_graph'`.
#' @param config A JSON formatted object containing the creation parameters for the output GraphML file. 
#'   Default is `{}`.
#' @param labels A vector of labels to assign to the name property of the vertex for display in cytoscape.
#'    If label value is NA, the mean mz value of the contracted vertex is assigned to the name property.
#'
#' @return No return value. Writes the generated similarity network as a GraphML file to
#'   the specified path.
#'
#' @details
#' This function processes a similarity matrix and corresponding spectral metadata to
#' create a graph object. It then exports the graph as a GraphML file, which is compatible
#' with popular network visualization tools like Cytoscape. The resulting network captures
#' the similarity between spectra based on fragment mass-to-charge ratios and MS1 annotations.
#'
#' @examples
#' # Example usage of outputGraphML function
#' outputGraphML(tbl = similarity_matrix, pmzs = precursor_mz, source = file_sources, 
#'               ioncount = ion_counts, annotation = annotations, 
#'               direction = directions, rt = retention_times, 
#'               contract = merge_indices, graph_name = "molecular_network")
#'
#' @import igraph
outputGraphML <- function(tbl, pmzs, source, ioncount, annotation, rt, contract, graph_name = 'new_graph', config = '{}', labels) {

  # Create a graph from a similarity matrix
  g <- igraph::graph.adjacency(
    1 - tbl,            # Convert similarity scores (1 = identical) to distances (0 = identical)
    mode     = "upper", # only uses top-right corner of table
    weighted = TRUE, 
    diag     = FALSE    # ignores the same-to-same values
  )

  # Assign vertex attributes
  V(g)$scans       <- V(g)$name # comes from the tbl row/col names
  V(g)$pmz         <- pmzs
  V(g)$source      <- source
  V(g)$ioncount    <- ioncount
  V(g)$rt          <- rt
  V(g)$avg_mz      <- pmzs  # placeholder for later averaging

  # Simplify the graph (remove loops and multiple edges)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  # Convert edge weights to absolute values
  E(g)$weight   <- abs(E(g)$weight)
  E(g)$ms1      <- ''

  # loop through any ms1 annotations and identify the vertices
  for (i in seq_along(annotation)) {
    transitions <- strsplit(annotation[i], ",")[[1]]

    for (trans in transitions) {
      side <- strsplit(trans, "->")[[1]]

      if (length(side) == 2) {
        vertex_info <- lapply(side, function(s) strsplit(s, ":")[[1]])

        if (all(sapply(vertex_info, length) == 2)) {
          vertices <- sapply(vertex_info, `[`, 1)
          ions     <- sapply(vertex_info, `[`, 2)

          # Assign the ms1 label to the matching edge
          E(g)[vertices[1] %--% vertices[2]]$ms1 <- paste0(ions[1], "->", ions[2])
        }
      }
    }
  }

  # Remove edges with original similarity = 0 (i.e., weight == 1 after conversion)
  g <- igraph::delete_edges(g, E(g)[which(E(g)$weight == 1)])

  # Contract vertices based on provided equivalence mapping
  g <- igraph::contract(g, contract, 
    vertex.attr.comb = 
    c(source   = unique, 
      avg_mz   = mean, 
      pmz      = toString, 
      rt       = toString, 
      scans    = toString, 
      ioncount = toString))

  # Simplify again after contraction
  g <- igraph::simplify(g, 
    remove.multiple = TRUE, 
    remove.loops    = TRUE, 
    edge.attr.comb  = list(weight = "min", ms1 = unique))

  # Set visual attributes (optional)
  V(g)$shape              <- "sphere"
  V(g)$color              <- "skyblue"
  V(g)$vertex.frame.color <- "white"

  # set to disply the mz values in cytoscape by default
  V(g)$name <- as.character(V(g)$avg_mz)

  # unless there is a designated label assigned to the scan
  hasLabel <- which(!is.na(labels))
  V(g)[contract[hasLabel]]$name <- labels[hasLabel]

  # Add graph-level attributes (key-value pairs)
  graph_attr(g, "generated_by") <- "NeutralLossGraph"
  graph_attr(g, "method") <- "graph.adjacency"
  graph_attr(g, "date") <- as.character(Sys.Date())
  graph_attr(g, "config") <- config

  # Export to GraphML
  igraph::write_graph(g, paste0(graph_name, ".graphml"), format = "graphml")
}
