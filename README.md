# Gene Expression Normalization Paper

Custom code for reproducing the results of the research paper by 
Roca, Gomes, Amorim & Scott-Fordsmand: "Variation-preserving normalization 
unveils blind spots in gene expression profiling".


## Instructions

To reproduce the results of the paper, please

1. Install/Update R and Bioconductor. Results were obtained with R version 
3.3.1 and Bioconductor version 3.3.

2. Install/Update required libraries: plotrix and RColorBrewer (R), and affy, 
drosgenome1.db, drosophila2.db, genefilter, and limma (Bioconductor).

3. Download the files of this repository into a local directory (or clone/fork 
this repository).

4. Download the raw microarray data of the *E. crypticus* dataset, and copy the 
`*.txt` files into the `rawdata` subdirectory. The dataset is available at 
NCBI/GEO, with GEO accessions 
[GSE69746](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69746), 
[GSE69792](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69792), 
[GSE69793](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69793) and
[GSE69794](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69794).

5. Download the raw microarray data of the Golden Spike dataset, and copy the 
`*.CEL` files into the `rawdata` subdirectory. The dataset is available at the 
online version of the paper by 
[Choe et al., *Genome Biology* 6, R16, 2005](http://dx.doi.org/10.1186/gb-2005-6-2-r16), 
in the additional data files 6 and 7. 

6. Download the raw microarray data of the Platinum Spike dataset, and copy the 
`*.cel` files into the `rawdata` subdirectory. The dataset is available at 
NCBI/GEO, with accession 
[GSE21344](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21344).

7. Execute the five R scripts `analyze_gene_expression*.r` **from** the local 
directory with the repository. For example, 

        $ cd ~/gene-expr-norm-paper
        $ nice R CMD BATCH --quiet --vanilla analyze_gene_expression.r
        $ nice R CMD BATCH --quiet --vanilla analyze_gene_expression_golden_spike.r
        $ nice R CMD BATCH --quiet --vanilla analyze_gene_expression_golden_spike_all_probe.r
        $ nice R CMD BATCH --quiet --vanilla analyze_gene_expression_platinum_spike.r
        $ nice R CMD BATCH --quiet --vanilla analyze_gene_expression_platinum_spike_all_probe.r


## Results

Execution of the scripts `analyze_gene_expression*.r` generates the following 
subdirectories:

* `bc_variation` between-condition variation (figures 2, S1)
* `boxplot` boxplots of expression levels (figure 1)
* `data` intermediate results
* `diff_expr` differential gene expression for the E. crypticus and synthetic 
datasets (figures 3-7, S2)
* `diff_expr_golden_spike` ROC curves for the Golden Spike dataset (figures 8, 
S4)
* `diff_expr_golden_spike_all_probe` ROC curves for the Golden Spike dataset 
without restricting probes (figure S3)
* `diff_expr_platinum_spike` ROC curves for the Platinum Spike dataset (figures 
9, S6)
* `diff_expr_platinum_spike_all_probe` ROC curves for the Platinum Spike 
dataset without restricting probes (figure S5)
* `pvalue` identification of no-variation genes (movies S4-S6)
* `stdvec` convergence of standard-vector normalization (movies S1-S3)
* `table` tables of experimental conditions for the E. crypticus dataset 
(tables S1, S2)

