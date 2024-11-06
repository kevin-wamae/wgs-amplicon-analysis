## workflow: wf-bcftools-illumina-amplicon

The wf-bcftools-illumina-amplicon directory contains scripts and resources for processing Illumina amplicon data using BCFtools. This directory is structured to facilitate the analysis of genomic data, including importing, processing, and summarizing variant data.

### Directory Structure

```
wf-bcftools-illumina-amplicon/
|-- input/
|   |-- study/
|       |-- extracted_variants/
|       |-- output/
|       |-- stats/
|-- scripts/
|   |-- bcftools-illumina.R
|-- R.project.Rproj
```

### Subdirectories and Files

- **`input/`**: Contains input data and directories for the analysis.
  - **`study/`**: Directory for storing study-specific files, such as sample names.
    - **`extracted_variants/`**: Directory for storing extracted variant files in TSV format, that have been generated using this Snakemake workflow -> [snakemake-illumina-bcftoolsvariant](https://github.com/kevin-wamae/snakemake-illumina-bcftoolsvariant)
    - **`output/`**: Directory for storing output files from the analysis.
    - **`stats/`**: Directory for storing statistics files related to read quality and mapping, also generated using the Snakemake workflow above.

- **`scripts/`**: Contains the main R script for processing the data.
  - **`main-bcftools.R`**: The main R script for processing Illumina amplicon data using BCFtools.

- **`R.project.Rproj`**: RStudio project file for this workflow.

### Dependencies

The scripts in this directory require R (version 4.0.0 or higher) and several R packages, including `tidyverse`, `data.table`, and `janitor`. Ensure you have these packages installed before running the scripts.
