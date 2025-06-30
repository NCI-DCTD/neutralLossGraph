# 🧪 msannotator

**Neutral Loss Graph** is an R package for generating annotated metabolomics data graphs using shared neutral losses. It simplifies the creation of GraphML outputs from MS features for downstream visualization and analysis.

---

## 📦 Installation

**Requires:** R version ≥ 4.3.0

You can install the package directly from GitHub using `devtools`:

```r
devtools::install_github("nci-dctd/neutralLossGraph")
```

## 🚀 Quick Start

Load the package and run a test graph generation:

```r
library(neutralLossGraph)
graphMinalemine()
```

This will generate a file in your current working directory:

- **File:** `package_test.graphml`

## 📁 Output

- **File:** `package_test.graphml`  
- **Location:** Current working directory  
- **Format:** [GraphML](https://en.wikipedia.org/wiki/GraphML)  
- **Description:**  
  A sample graph generated from internal test data representing metabolomic feature annotations.  
  This file can be visualized using tools like [Cytoscape](https://cytoscape.org/), Gephi, or other GraphML-compatible software.
