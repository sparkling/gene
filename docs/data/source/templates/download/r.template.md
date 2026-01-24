## R Package

### Installation

```r
# From CRAN
install.packages("{{PACKAGE_NAME}}")

# From Bioconductor
BiocManager::install("{{PACKAGE_NAME}}")

# From GitHub
devtools::install_github("{{GITHUB_REPO}}")
```

### Quick Start

```r
library({{package_name}})

# Load data
data <- {{function_name}}({{PARAMS}})

# Process
result <- {{process_function}}(data)
```

### Dependencies

| Package | Source | Purpose |
|---------|--------|---------|
| {{DEP_1}} | {{CRAN/Bioc}} | {{PURPOSE}} |
| {{DEP_2}} | {{CRAN/Bioc}} | {{PURPOSE}} |

### Common Operations

#### Operation 1: {{OPERATION_NAME}}

```r
{{CODE_EXAMPLE_1}}
```

#### Operation 2: {{OPERATION_NAME}}

```r
{{CODE_EXAMPLE_2}}
```

### Vignettes

```r
browseVignettes("{{PACKAGE_NAME}}")
```

### Documentation

{{PACKAGE_DOCS_URL}}
