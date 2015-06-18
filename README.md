# Gene Expression Normalization Paper

Custom code for reproducing the results of the research paper by 
Roca, Gomes, Amorim & Scott-Fordsmand: "A novel normalization approach 
unveils blind spots in gene expression profiling". 


## Instructions

To reproduce the results of the paper, please

1. Install/Update R and Bioconductor. Results were obtained with R version 
3.2.0 and Bioconductor version 3.1.

2. Install required libraries: plotrix and RColorBrewer (R), genefilter and 
limma (Bioconductor).

3. Download the files in this repository into a local directory (or fork/clone 
this repository).

4. Download the raw microarray data from NCBI/GEO into a **single** local 
directory. This directory can (and should) be different from the script 
directory (previous step). The GEO accessions are 
[GSE69746](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69746), 
[GSE69792](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69792), 
[GSE69793](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69793) and
[GSE69794](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69794).

5. Edit the line 45 of file `analyze_gene_expression.r` so that the variable 
`ma.data.dir` contains the name of the raw data directory (previous step). For 
example, change

        ma.data.dir <- ""   # raw data directory
into

        ma.data.dir <- "~/gene-expr-norm-paper/rawdata"

6. Execute the R script `analyze_gene_expression.r` from its directory. For 
example, 

        $ cd ~/gene-expr-norm-paper/code
        $ nice R CMD BATCH --quiet --vanilla analyze_gene_expression.r &


## Results

Execution of the script `analyze_gene_expression.r` generates the following 
subdirectories:
- `bc_variation`  between-condition variation (figure 2)
- `boxplot`       boxplots of expression levels (figure 1)
- `data`          intermediate results
- `diff_expr`     differential expression (figure 3, extended data figures 1-2)
- `pvalue`        identification of non-differentially expressed genes (videos 
4-5)
- `stdvec`        convergence of standard-vector normalization (videos 1-3)
- `table`         experimental condition tables (extended data tables 1-2)


