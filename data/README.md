## Data

- **`01_numbering_reference_strain/`**  
  Sequence of vaccine strain used as reference for defining the location of epitope file in `R/config.R`.

- **`02_vaccine_strain/`**  
  Contain all the vaccine strains you need to calculate antigenicity. It is already one year. You can update the vaccine strain if it is announced and reported in the GISAID database.

- **`03_tmp/`**  
  This is the demo sample. It will have FASTQ for demo data and, importantly, the HA structure in PDB for Flu A H1N1 (I think). You may need to include other subtypes in the dashboard too.

- **`04_circulating_strain/`**  
  Provide all circulating strains since 20XX something to 2025. This will be used in the dashboard to compare with input FASTQ and show the user how much the targeted strain escapes the vaccine compared with circulating strains.
