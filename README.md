## Project Overview

This repository contains scripts and resources for genomic data analysis, specifically focusing on the study of malaria parasites.

- The project includes various workflows for processing and analyzing genomic data from different sequencing platforms, such as Illumina and Nanopore.
- The primary goal is to aggregate and analyze allele frequency data, compute haplotype frequencies, and perform quality control on the extracted data.

## Directory Structure

The repository is organized into the following directories:

```
.
|-- readme.md
|-- .gitignore
|-- resources-genome
|-- wf-bcftools-illumina-amplicon
|-- wf-nanorave-nanopore-amplicon
`-- wf-seekdeep-amplicon
```

- [**resources-genome/**](resources-genome/README.md): Contains genomic resources, including BED and FASTA files.
- [**wf-bcftools-illumina-amplicon/**](wf-bcftools-illumina-amplicon/README.md): Contains scripts and resources for processing Illumina amplicon data using BCFtools.
- [**wf-nanorave-nanopore-amplicon/**](wf-nanorave-nanopore-amplicon/README.md): Contains scripts and resources for processing Nanopore amplicon data using NanoRave.
- [**wf-seekdeep-amplicon/**](wf-seekdeep-amplicon/README.md): Contains scripts and resources for processing amplicon data (Illumina and Nanopore) using SeekDeep.

## Getting Started

To get started with this project, clone the repository and navigate to the relevant workflow directory based on your sequencing platform. Each workflow directory contains detailed scripts and instructions for processing and analyzing the genomic data.

## Dependencies

The scripts in this repository require R (version 4.0.0 or higher) and several R packages, including `tidyverse`, `data.table`, and `janitor`. Ensure you have these packages installed before running the scripts.

## Usage

1. Clone the repository:

    ```sh
    git clone https://github.com/kkariuki/wgs-amplicon-analysis.git
    ```

2. Navigate to the relevant workflow directory and read the instructions:

3. Follow the instructions in the scripts to process and analyze your genomic data.

## **Contributing and Issues**

Contributions to this project are welcome. Please fork the repository and submit a pull request with your changes.

Additionally, report any issues or bugs by openning an issue [here](https://github.com/kkariuki/wgs-amplicon-analysis/issues) or contact me via email at **wamaekevin[at]gmail[dot]com**

## License

This project is licensed under the MIT License. See the LICENSE file for details.
