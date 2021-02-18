
## Master R script for:

### Changes in the Blood Transcriptome Following Treatment Reflect Latent Tuberculosis Heterogeneity

*Julie G Burel, Akul Singhania, Paige Dubelko, Julius Muller, Rachel Tanner, Eneida Parizotto, Martin Dedicoat, Thomas E. Fletcher, James Dunbar, Adam F. Cunningham, Cecilia S. Lindestam Arlehamn, Donald Catanzaro, Antonino Catanzaro, Timothy Rodwell, Helen McShane, Matthew K. O’Shea and Bjoern Peters*

### Burel_et_al_R_scripts.R
(author - Akul Singhania)
The following steps are outlined in this script:
1) Read in raw Illumina IDAT files
2) Process raw data into a normalized matrix
3) Quality check
4) Filtering the data
5) Annotation of Illumina ID probes to Gene symbols
6) Principal Component Analysis (PCA)
7) Batch correction
8) Differential gene expression analysis (DEG)

### compileAnnotations.R 
(author - Paige Dubelko)
First script for getting Illumina Annotations 

Annotations for Illumina microarrays were collected from three sources:
1) BioMart - collected using an R package
2) GEO - Illumina specific annotations downloaded from NCBI GEO database
3) Illumina - Annotation downloaded from the Illumina website.

This script then merges together the different annotation sources.

### geneAnnotation.py 
(author - Paige Dubelko)
Second script for getting Illumina Annotations
