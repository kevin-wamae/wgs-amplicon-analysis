## resources-genome/

The resources-genome directory contains genomic resources essential for the analysis of malaria parasites. This directory includes various file formats such as BED and FASTA files, which are used for storing genomic sequences and annotations.

### Directory Structure

```
resources-genome/
|-- bed/
|   |-- Pf3D7_v51.bed
|-- fasta-protein/
|   |-- PfCRT.txt
|   |-- PfDHFR.txt
|   |-- PfDHPS.txt
|   |-- PfK13.txt
|   |-- PfMDR1.txt
```

### Subdirectories and Files

- **`bed/`**: This subdirectory contains BED files, which are used to store genomic features such as genes, exons, and other annotations.
  - **`Pf3D7_v51.bed`**: A BED file containing annotations for the Plasmodium falciparum genome, version 51. This file includes information about various protein-coding genes and their descriptions and can be obtained from the [PlasmoDB](https://plasmodb.org) database.

- **`fasta-protein/`**: This subdirectory contains FASTA files, which are used to store protein sequences.
  - **`PfCRT.txt`**: A FASTA file containing a protein sequence.

### Example BED File Content

Here is an example excerpt from the

Pf3D7_v51.bed

 file:

```bed
Pf3D7_05_v3 1013253 1015395 PF3D7_0524400 . - VEuPathDB protein_coding_gene . ID=PF3D7_0524400;description=ribosome-interacting GTPase 1%2C putative
Pf3D7_05_v3 1015607 1017341 PF3D7_0524500 . + VEuPathDB protein_coding_gene . ID=PF3D7_0524500;description=conserved Plasmodium protein%2C unknown function
Pf3D7_05_v3 1017357 1019072 PF3D7_0524600 . + VEuPathDB protein_coding_gene . ID=PF3D7_0524600;description=50S ribosomal protein L12%2C apicoplast%2C putative
```

### Example FASTA File Content

Here is an example excerpt from the PfCRT FASTA file (*note that this is not the standard FASTA format which starts with a `>` character*):

```fasta
MKFASKKNNQKNSSKNDERYRELDNLVQEGNGSRLGGGSCLGKCAHVFKLIFKEIKDNIFIYILSIIYLSVCVMNKIFAKRTLNKIGNYSFVTSETHNFICMIMFFIVYSLFGNKKGNSKERHRSFNLQFFAISMLDACSVILAFIGLTRTTGNIQSFVLQLSIPINMFFCFLILRYRYHLYNYLGAVIIVVTIALVEMKLSFETQEENSIIFNLVLISALIPVCFSNMTREIVFKKYKIDILRLNAMVSFFQLFTSCLILPVYTLPFLKQLHLPYNEIWTNIKNGFACLFLGRNTVVENCGLGMAKLCDDCDGAWKTFALFSFFNICDNLITSYIIDKFSTMTYTIVSCIQGPAIAIAYYFKFLAGDVVREPRLLDFVTLFGYLFGSIIYRVGNIILERKKMRNEENEDSEGELTNVDSIITQ
```

These files are crucial for genomic data analysis, providing the necessary annotations and sequences for various computational workflows.
