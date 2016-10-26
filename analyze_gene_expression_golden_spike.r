# Copyright (c) 2016, Universitat Rovira i Virgili (Spain), Aarhus University 
# (Denmark) and University of Aveiro (Portugal)
# 
# Written by Carlos P. Roca
# as Research funded by the European Union
# for the research paper by Roca, Gomes, Amorim & Scott-Fordsmand: "Variation-
# preserving normalization unveils blind spots in gene expression profiling".
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# Analyzes gene expression in golden spike dataset


library( affy )
library( genefilter )
library( limma )
library( RColorBrewer )

source( "normalize_median_condec.r" )
source( "normalize_stdvec_condec.r" )
source( "select_nodup.r" )
source( "watson_u2.r" )


# set directories

ma.data.dir <- ""   # raw data directory

stopifnot( ma.data.dir != "" )

gene.diff.expr.dir <- "gene_diff_expr_golden_spike"


# define experimental conditions and samples

ma.treatment <- c( "spike" )
ma.control <- c( "constant" )
ma.condition <- c( ma.treatment, unique( ma.control ) )

ma.condition.sample.num <- c( 3, 3 )
ma.sample.condition <- rep( ma.condition, ma.condition.sample.num )
ma.sample.postfix <- unlist( lapply( ma.condition.sample.num, 
    function( mcsn ) 1:mcsn ) )
ma.sample <- paste0( ma.sample.condition, '.', ma.sample.postfix )

ma.sample.num <- length( ma.sample )
ma.condition.num <- length( ma.condition )
ma.treatment.num <- length( ma.treatment )


# read platform annotation and identify probes

ma.species <- "drosgenome1"
ma.db.name <- paste0( ma.species, ".db" )

library( ma.db.name, character.only=TRUE )
ma.db <- get( ma.db.name )

ma.probe.platform <- sort( keys( ma.db, "PROBEID" ) )

ma.annot <- select.nodup( ma.db, ma.probe.platform, 
    c( "SYMBOL", "GENENAME", "FLYBASE" ), "PROBEID" )

names( ma.annot ) <- sub( "PROBEID", "probe.id", names( ma.annot ), fixed=TRUE )


# read probe set fold change from experiment design

exp.design.file.name <- "golden_spike_probe_set_fold_change.csv"

exp.design <- read.csv( exp.design.file.name, stringsAsFactors=FALSE, 
    comment.char="#" )

names( exp.design ) <- sub( "probe.set.id", "probe", names( exp.design ), 
    fixed=TRUE )
names( exp.design ) <- sub( "fold.change", "fc", names( exp.design ), 
    fixed=TRUE )

stopifnot( exp.design$probe == ma.probe.platform )


# read raw expression data

ma.data.def.file.name <- "golden_spike_data_def.txt"

ma.data.def <- read.delim( ma.data.def.file.name, stringsAsFactors=FALSE )

ma.sample.file.name.read <- ma.data.def$FileName
ma.sample.read <- ma.data.def$Sample
ma.sample.condition.read <- ma.data.def$Condition

ma.raw.data <- ReadAffy( filenames = 
        file.path( ma.data.dir, ma.data.def$FileName ) )


# verify probes of raw data

stopifnot( sort( featureNames( ma.raw.data ) ) == ma.probe.platform )


# verify samples of raw data

stopifnot( sampleNames( ma.raw.data ) == ma.sample.file.name.read )

stopifnot( sort( ma.sample ) == sort( ma.sample.read ) )

stopifnot( ma.sample.condition[ order( ma.sample ) ] == 
    ma.sample.condition.read[ order( ma.sample.read ) ] )


# get background-corrected data

ma.expr.data <- log2( exprs( mas5( ma.raw.data, normalize=FALSE ) ) )


# remove expression data below noise threshold

ma.mas5.call <- exprs( mas5calls( ma.raw.data ) )
ma.expr.data.missing.bol <- ma.mas5.call != "P"

stopifnot( rownames( ma.expr.data ) == rownames( ma.expr.data.missing.bol ) )
stopifnot( colnames( ma.expr.data ) == colnames( ma.expr.data.missing.bol ) )

ma.expr.data[ ma.expr.data.missing.bol ] <- NA


# remove control probes and probes with missing data
# reorder probes in expression data

ma.probe.control.bol <- grepl( "^AFFX-", ma.probe.platform )
ma.probe.missing.data.bol <- rowSums( is.na( ma.expr.data ) ) > 0

ma.probe <- ma.probe.platform[ ! ma.probe.control.bol &
    ! ma.probe.missing.data.bol ]
ma.probe.num <- length( ma.probe )

ma.expr.data <- ma.expr.data[ ma.probe, ]


# identify probes from experiment design

ma.probe.known.pos <- intersect( ma.probe, 
    exp.design$probe[ exp.design$fc > 1 ] )
ma.probe.known.neg <- intersect( ma.probe, 
    exp.design$probe[ exp.design$fc == 1 ] )
ma.probe.known <- intersect( ma.probe, 
    exp.design$probe[ exp.design$fc >= 1 ] )
ma.probe.unknown <- intersect( ma.probe, 
    exp.design$probe[ exp.design$fc < 1 ] )


# reorder and rename samples in expression data

ma.expr.data <- ma.expr.data[ , match( ma.sample, ma.sample.read ) ]

colnames( ma.expr.data ) <- ma.sample


# identify samples and conditions for analysis

set.treatment <- ma.treatment
set.control <- ma.control
set.condition <- ma.condition

set.sample <- ma.sample
set.sample.condition <- ma.sample.condition

set.treatment.num <- length( set.treatment )
set.condition.num <- length( set.condition )
set.sample.num <- length( set.sample )


# define comparisons between experimental conditions

treatment.control.comparison <- cbind( set.treatment, set.control )

set.comparison <- rbind( treatment.control.comparison )

colnames( set.comparison ) <- NULL

set.comparison.name <- apply( set.comparison, 1, paste0, collapse=".vs." )
set.comparison.num <- length( set.comparison.name )


# select expression data for analysis

set.expr.data <- ma.expr.data[ , ma.sample.condition %in% set.condition ]


# normalize expression data

random.seed <- 20160801
set.seed( random.seed )

norm.method <- c( "no", "median", "quantile", "theoretical.loess", "loess", 
    "theoretical.median.condec", "median.condec", "theoretical.stdvec.condec", 
    "stdvec.condec" )

norm.method.pch <- c( 8, 1, 2, 23, 23, 21, 21, 24, 24 )
norm.method.lty <- c( 1, 1, 1, 2, 1, 2, 1, 2, 1 )
norm.method.col <- c( "orange3", "black", "red3", "magenta3", "magenta3", 
    "green3", "green3", "blue3", "blue3" )

bg.col.offset <- c( 0.8, 0.8, 0.8, 0 )
norm.method.bg.col <- c( "orange3", "black", "red3", 
    adjustcolor( "magenta3", offset=bg.col.offset ), "magenta3", 
    adjustcolor( "green3", offset=bg.col.offset ), "green3", 
    adjustcolor( "blue3", offset=bg.col.offset ), "blue3" )

names( norm.method.pch ) <- norm.method
names( norm.method.lty ) <- norm.method
names( norm.method.col ) <- norm.method
names( norm.method.bg.col ) <- norm.method

for ( nmeth in norm.method )
{
    if ( nmeth == "no" )
    {
        no.norm.data <- set.expr.data
    }
    else if ( nmeth == "median" )
    {
        median.offset <- apply( set.expr.data, 2, median, na.rm=TRUE )
        median.offset <- median.offset - mean( median.offset )
        median.norm.data <- sweep( set.expr.data, 2, median.offset )
    }
    else if ( nmeth == "quantile" )
    {
        quantile.norm.data <- 
            normalizeBetweenArrays( set.expr.data, method="quantile" )
    }
    else if ( nmeth == "theoretical.loess" )
    {
        theoretical.loess.norm.data <- 
            normalize.loess( set.expr.data, 
                subset = match( ma.probe.known.neg, ma.probe ) )
    }
    else if ( nmeth == "loess" )
    {
        loess.norm.data <- normalize.loess( set.expr.data )
    }
    else if ( nmeth == "theoretical.median.condec" )
    {
        theoretical.median.condec.norm.result <- 
            normalize.median.condec( set.expr.data, 
            set.sample.condition, ma.probe.known.neg, search.h0.probe=FALSE, 
            verbose=TRUE )
        theoretical.median.condec.norm.data <- 
            theoretical.median.condec.norm.result$data
    }
    else if ( nmeth == "median.condec" )
    {
        median.condec.norm.result <- normalize.median.condec( set.expr.data, 
            set.sample.condition, convergence.threshold = c( 0.001, 0.3 ), 
            verbose=TRUE )
        median.condec.norm.data <- median.condec.norm.result$data
    }
    else if ( nmeth == "theoretical.stdvec.condec" )
    {
        theoretical.stdvec.condec.norm.result <- 
            normalize.stdvec.condec( set.expr.data, 
            set.sample.condition, ma.probe.known.neg, search.h0.probe=FALSE, 
            verbose=TRUE )
        theoretical.stdvec.condec.norm.data <- 
            theoretical.stdvec.condec.norm.result$data
    }
    else if ( nmeth == "stdvec.condec" )
    {
        stdvec.condec.norm.result <- normalize.stdvec.condec( set.expr.data, 
            set.sample.condition, verbose=TRUE )
        stdvec.condec.norm.data <- stdvec.condec.norm.result$data
    }
    else
        stop( "wrong normalization method ", nmeth )
}


# build model matrix for analysis of differential gene expression with limma

condition.factor <- factor( set.sample.condition, levels=set.condition )

model.design <- model.matrix( ~0 + condition.factor )
rownames( model.design ) <- set.sample
colnames( model.design ) <- set.condition

diff.expr.contrast <- apply( set.comparison, 1, function( sc ) {
    de.cont <- rep( 0, set.condition.num )
    de.cont[ set.condition == sc[ 1 ] ] <- 1
    de.cont[ set.condition == sc[ 2 ] ] <- -1
    de.cont
} )

rownames( diff.expr.contrast ) <- set.condition
colnames( diff.expr.contrast ) <- set.comparison.name


# analyze differential gene expression by iterating first over normalization 
# methods and then over differential gene expression methods

gene.diff.expr.method <- c( "ttest", "limma" )

gene.test.alpha <- 0.05

set.comparison.gene.test <- as.vector( sapply( set.comparison.name, paste0, 
    ".", c( "fc", "pv", "fdr" ) ) )

for ( nmeth in norm.method )
{
    norm.data <- get( paste0( nmeth, ".norm.data" ) )
    
    stopifnot( rownames( norm.data ) == ma.probe )
    stopifnot( colnames( norm.data ) == set.sample )
    
    for ( gdemeth in gene.diff.expr.method )
    {
        if ( gdemeth == "ttest" )
        {
            gene.diff.expr <-  apply( set.comparison, 1, function( sc ) {
                tc.col <- set.sample.condition == sc[ 1 ] | 
                    set.sample.condition == sc[ 2 ]
                tc.factor <- factor( ( 1*( set.sample.condition == sc[ 1 ] ) + 
                        2*( set.sample.condition == sc[ 2 ] ) )[ tc.col ] )
                rtt <- rowttests( norm.data[ , tc.col ], tc.factor )
                rtt$fdr <- p.adjust( rtt$p.value, method="fdr" )
                c( rtt$dm, rtt$p.value, rtt$fdr )
            } )
        }
        else if ( gdemeth == "limma" )
        {
            model.fit <- lmFit( norm.data, model.design )
            diff.expr.fit <- contrasts.fit( model.fit, diff.expr.contrast )
            diff.expr.efit <- eBayes( diff.expr.fit )
            gene.diff.expr <- sapply( set.comparison.name, function( scn ) {
                de.table <- topTable( diff.expr.efit, coef=scn, number=Inf, 
                    sort.by="none", adjust.method="fdr" )
                c( de.table$logFC, de.table$P.Value, de.table$adj.P.Val )
            } )
        }
        else
            stop( "wrong gene differential expression method ", gdemeth )
        
        dim( gene.diff.expr ) <- c( ma.probe.num, 
            length( set.comparison.gene.test ) )
        rownames( gene.diff.expr ) <- ma.probe
        colnames( gene.diff.expr ) <- set.comparison.gene.test
        
        gene.diff.expr <- data.frame( probe=ma.probe, 
            gene.diff.expr[ ma.probe, ], stringsAsFactors=FALSE )
        
        assign( paste0( nmeth, ".norm.", gdemeth, ".gene.diff.expr" ), 
            gene.diff.expr )
    }
}


# plot roc curves

plot.norm.method <- c( "median", "quantile", "theoretical.loess", "loess", 
    "theoretical.median.condec", "median.condec", "theoretical.stdvec.condec", 
    "stdvec.condec" )

plot.fdr.point <- c( 0.01, 0.05 )

if ( ! file.exists( gene.diff.expr.dir ) )
    dir.create( gene.diff.expr.dir )

for ( gdemeth in gene.diff.expr.method )
{
    rp.file.name <- paste0( gdemeth, ".roc.pdf" )
    cairo_pdf( filename = file.path( gene.diff.expr.dir, rp.file.name ), 
        width=3.25, height=2.60 )
    
    par( mar = c( 1.90, 1.70, 0.15, 0.15 ), mgp = c( 0.90, 0, 0 ), tcl=-0.10 )
    
    plot( 1, type="n", cex.lab=0.67, cex.axis=0.67, 
        xlim = c( 0, 1 ), ylim = c( 0, 1 ), 
        xlab="False positive rate", 
        ylab="True positive rate" )
    
    for ( nmeth in plot.norm.method )
    {
        gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
            ".gene.diff.expr" ) )
        
        gene.diff.expr.probe <- gene.diff.expr$probe
        gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
        
        gene.diff.expr.probe <- gene.diff.expr.probe[ 
            order( gene.diff.expr.fdr ) ]
        gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
        
        gene.diff.expr.fprtpr <- sapply( 1 : length( gene.diff.expr.probe ), 
            function( i ) {
                gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                c( mean( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                    mean( ma.probe.known.pos %in% gene.diff.expr.pos ) )
            } )
        
        lines( gene.diff.expr.fprtpr[ 1, ], gene.diff.expr.fprtpr[ 2, ], lwd=1, 
            lty = norm.method.lty[ nmeth ], col = norm.method.col[ nmeth ] )
        
        fdr.point <- sapply( plot.fdr.point, function( pfp ) 
            sum( gene.diff.expr.fdr < pfp ) )
        
        points( gene.diff.expr.fprtpr[ 1, fdr.point ], 
            gene.diff.expr.fprtpr[ 2, fdr.point ], 
            pch = norm.method.pch[ nmeth ], col = norm.method.col[ nmeth ], 
            bg = norm.method.bg.col[ nmeth ] )
    }
    
    dev.off()
}


# plot roc curves near origin

plot.fdr.point.ref <- c( plot.fdr.point, 0.1, 0.2, 0.5 )

for ( gdemeth in gene.diff.expr.method )
{
    rpo.file.name <- paste0( gdemeth, ".roc.origin.pdf" )
    cairo_pdf( filename = file.path( gene.diff.expr.dir, rpo.file.name ), 
        width=3.25, height=2.60 )
    
    par( mar = c( 1.90, 1.70, 0.15, 0.15 ), mgp = c( 0.90, 0, 0 ), tcl=-0.10 )
    
    plot( 1, type="n", cex.lab=0.67, cex.axis=0.67, 
        xlim = c( 0, 870 ), ylim = c( 0, 1150 ), 
        xlab="# False positives", 
        ylab="# True positives" )
    
    for ( pfpr in plot.fdr.point.ref )
        abline( a=0, b = (1-pfpr)/pfpr, lty=2, lwd=1 )
    
    for ( nmeth in plot.norm.method )
    {
        gene.diff.expr <- get( paste0( nmeth, ".norm.", gdemeth, 
            ".gene.diff.expr" ) )
        
        gene.diff.expr.probe <- gene.diff.expr$probe
        gene.diff.expr.fdr <- gene.diff.expr$spike.vs.constant.fdr
        
        gene.diff.expr.probe <- gene.diff.expr.probe[ 
            order( gene.diff.expr.fdr ) ]
        gene.diff.expr.fdr <- sort( gene.diff.expr.fdr )
        
        gene.diff.expr.fptp <- sapply( 1 : length( gene.diff.expr.probe ), 
            function( i ) {
                gene.diff.expr.pos <- gene.diff.expr.probe[ 1 : i ]
                c( sum( ma.probe.known.neg %in% gene.diff.expr.pos ), 
                    sum( ma.probe.known.pos %in% gene.diff.expr.pos ) )
            } )
        
        lines( gene.diff.expr.fptp[ 1, ], gene.diff.expr.fptp[ 2, ], lwd=1, 
            lty = norm.method.lty[ nmeth ], col = norm.method.col[ nmeth ] )
        
        fdr.point <- sapply( plot.fdr.point, function( pfp ) 
            sum( gene.diff.expr.fdr < pfp ) )
        
        points( gene.diff.expr.fptp[ 1, fdr.point ], 
            gene.diff.expr.fptp[ 2, fdr.point ], 
            pch = norm.method.pch[ nmeth ], col = norm.method.col[ nmeth ],
            bg = norm.method.bg.col[ nmeth ] )
    }
    
    dev.off()
}


# print probe numbers

cat( sprintf( "number of platform probes: %d\n", 
    sum( ! ma.probe.control.bol ) ) )
cat( sprintf( "number of spike-ins: %d\n", 
    sum( exp.design$fc != -1 & ! ma.probe.control.bol ) ) )
cat( sprintf( "number of positive spike-ins: %d\n", 
    sum( exp.design$fc > 1 & ! ma.probe.control.bol ) ) )
cat( sprintf( "number of negative spike-ins: %d\n", 
    sum( exp.design$fc == 1 & ! ma.probe.control.bol ) ) )

cat( sprintf( "number of probes: %d\n", ma.probe.num ) )
cat( sprintf( "number of unknown probes: %d\n", length( ma.probe.unknown ) ) )
cat( sprintf( "number of known probes: %d\n", length( ma.probe.known ) ) )
cat( sprintf( "number of known positives: %d\n", 
    length( ma.probe.known.pos ) ) )
cat( sprintf( "number of known negatives: %d\n", 
    length( ma.probe.known.neg ) ) )

h0.probe.norm.method <- c( "median.condec", "stdvec.condec" )

for ( hpnmeth in h0.probe.norm.method )
{
    norm.result <- get( paste0( hpnmeth, ".norm.result" ) )
    norm.h0.probe <- norm.result$h0.probe
    
    cat( sprintf( "%s - h0 probes - number: %d\n", hpnmeth, 
        length( norm.h0.probe ) ) )
    cat( sprintf( "%s - h0 probes - fraction of known probes: %g\n", hpnmeth, 
        mean( norm.h0.probe %in% ma.probe.known ) ) )
    cat( sprintf( 
        "%s - h0 probes - fraction of known negatives among known probes: %g\n", 
        hpnmeth, sum( norm.h0.probe %in% ma.probe.known.neg ) / 
            sum( norm.h0.probe %in% ma.probe.known ) ) )
}

