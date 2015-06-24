# Copyright (c) 2015, Universitat Rovira i Virgili (Spain), University of 
# Aveiro (Portugal) & Aarhus University (Denmark)
# 
# Written by Carlos P. Roca
# as Research funded by the European Union
# for the research paper by Roca, Gomes, Amorim & Scott-Fordsmand: "A novel 
# normalization approach unveils blind spots in gene expression profiling"
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


# Analyzes gene expression variation using several normalization algorithms


library( genefilter )
library( limma )
library( RColorBrewer )


source( "normalize_stdvec_condec.r" )
source( "normalize_stdvec_condec_probe.r" )

source( "normalize_median_condec.r" )
source( "normalize_median_condec_probe.r" )


# set directories

ma.data.dir <- ""   # raw data directory

stopifnot( ma.data.dir != "" )

data.dir <- "data"

table.dir <- "table"
boxplot.dir <- "boxplot"
p.value.dir <- "pvalue"
stdvec.dir <- "stdvec"
bc.variation.dir <- "bc_variation"
diff.expr.dir <- "diff_expr"


# read experiment data

targets.file.name <- "ecrypticus_targets.txt"

targets <- readTargets( file=targets.file.name )

ma.sample.read <- targets$Label
ma.sample.condition.read <- targets$condition


# define experimental conditions

ag.treatment <- c( 
    "Ag.AgNO3.EC20.3d", "Ag.AgNO3.EC20.7d", 
    "Ag.AgNO3.EC50.3d", "Ag.AgNO3.EC50.7d", 
    "Ag.Coated.EC20.3d", "Ag.Coated.EC20.7d", 
    "Ag.Coated.EC50.3d", "Ag.Coated.EC50.7d", 
    "Ag.NC.EC20.3d", "Ag.NC.EC20.7d", 
    "Ag.NC.EC50.3d", "Ag.NC.EC50.7d", 
    "Ag.NM300K.EC20.3d", "Ag.NM300K.EC20.7d", 
    "Ag.NM300K.EC50.3d", "Ag.NM300K.EC50.7d"
)

ag.control <- c( 
    "Ag.CT.3d", "Ag.CT.7d", 
    "Ag.CT.3d", "Ag.CT.7d", 
    "Ag.CT.3d", "Ag.CT.7d", 
    "Ag.CT.3d", "Ag.CT.7d", 
    "Ag.CT.3d", "Ag.CT.7d", 
    "Ag.CT.3d", "Ag.CT.7d", 
    "Ag.CTD.3d", "Ag.CTD.7d", 
    "Ag.CTD.3d", "Ag.CTD.7d"
)
    
cu.treatment <- c( 
    "Cu.CuNO3.EC20.3d", "Cu.CuNO3.EC20.7d", 
    "Cu.CuNO3.EC50.3d", "Cu.CuNO3.EC50.7d", 
    "Cu.Cu.NPs.EC20.3d", "Cu.Cu.NPs.EC20.7d", 
    "Cu.Cu.NPs.EC50.3d", "Cu.Cu.NPs.EC50.7d", 
    "Cu.Cu.Nwires.EC20.3d", "Cu.Cu.Nwires.EC20.7d", 
    "Cu.Cu.Nwires.EC50.3d", "Cu.Cu.Nwires.EC50.7d", 
    "Cu.Cu.field.EC20.3d", "Cu.Cu.field.EC20.7d", 
    "Cu.Cu.field.EC50.3d", "Cu.Cu.field.EC50.7d"
)

cu.control <- c( 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d", 
    "Cu.CT.3d", "Cu.CT.7d" 
)

ni.treatment <- c( 
    "Ni.NiNO3.EC20.3d", "Ni.NiNO3.EC20.7d", 
    "Ni.NiNO3.EC50.3d", "Ni.NiNO3.EC50.7d", 
    "Ni.Ni.NPs.EC20.3d", "Ni.Ni.NPs.EC20.7d", 
    "Ni.Ni.NPs.EC50.3d", "Ni.Ni.NPs.EC50.7d"
)

ni.control <- c( 
    "Ni.CT.3d", "Ni.CT.7d", 
    "Ni.CT.3d", "Ni.CT.7d", 
    "Ni.CT.3d", "Ni.CT.7d", 
    "Ni.CT.3d", "Ni.CT.7d"
)

uv.treatment <- c( "Uv.UV.D1.5d", "Uv.UV.D2.5d" )

uv.control <- c( "Uv.CT.5d", "Uv.CT.5d" )

ag.condition <- c( ag.treatment, unique( ag.control ) )
cu.condition <- c( cu.treatment, unique( cu.control ) )
ni.condition <- c( ni.treatment, unique( ni.control ) )
uv.condition <- c( uv.treatment, unique( uv.control ) )

ma.treatment <- c( ag.treatment, cu.treatment, ni.treatment, uv.treatment )
ma.control <- c( ag.control, cu.control, ni.control, uv.control )
ma.condition <- c( ag.condition, cu.condition, ni.condition, uv.condition )

ma.sample <- as.vector( sapply( ma.condition, paste0, ".", 1:3 ) )
ma.sample.condition <- rep( ma.condition, each=3 )

# check matching with targets file
stopifnot( sort( ma.sample ) == sort( ma.sample.read ) )
stopifnot( sort( ma.sample.condition ) == sort( ma.sample.condition.read ) )

ma.sample.num <- length( ma.sample )
ma.condition.num <- length( ma.condition )
ma.treatment.num <- length( ma.treatment )

# write table of experimental conditions
if ( ! file.exists( table.dir ) )
    dir.create( table.dir )

write.csv( ma.condition, file.path( table.dir, "experimental_condition.csv" ) )


# define colors for experimental conditions

ma.element <- c( "ag", "cu", "ni", "uv" )

for ( mel in ma.element )
{
    mel.condition <- get( paste0( mel, ".condition" ) )
    mel.condition.num <- length( mel.condition )
    
    mel.color.brewer <- brewer.pal( 8, "Set1" )
    mel.color <- switch( mel, 
        ag = c( adjustcolor( mel.color.brewer[ 2 ], blue.f=1.3 ),   # blue
            mel.color.brewer[ 6 ] ),  # yellow
        cu = c( adjustcolor( mel.color.brewer[ 1 ], red.f=1.1 ),    # red
            rgb( 0.6, 1, 1 ) ),       # cyan
        ni = c( adjustcolor( mel.color.brewer[ 3 ], green.f=1.2 ),  # green
            mel.color.brewer[ 5 ] ),  # orange
        uv = c( mel.color.brewer[ 4 ], mel.color.brewer[ 8 ] )      # purple
    )
    
    mel.condition.color <- colorRampPalette( mel.color )( mel.condition.num )
    
    # reorder colors to maximize contrast between treatments of the same element
    mel.condition.color <- mel.condition.color[ sapply( 1:mel.condition.num, 
        function( i ) ifelse( i%%2 == 1, (i+1)/2, 
            i/2 + (mel.condition.num-1)%/%2 + 1 )
    ) ]
    
    names( mel.condition.color ) <- mel.condition
    
    assign( paste0( mel, ".condition.color" ), mel.condition.color )
}

ma.condition.color <- c( ag.condition.color, cu.condition.color, 
    ni.condition.color, uv.condition.color )

ma.sample.color <- ma.condition.color[ ma.sample.condition ]
names( ma.sample.color ) <- ma.sample

ma.treatment.color <- ma.condition.color[ ma.treatment ]


# read microarray data

current.dir <- setwd( ma.data.dir )

ma.data <- read.maimages( targets, source="agilent.median", green.only=TRUE, 
    other.columns = c( "gProcessedSignal", "gIsPosAndSignif" ) )

setwd( current.dir )


# remove levels according to agilent flag IsPosAndSignif

raw.expr.data <- ma.data$other$gProcessedSignal
row.names( raw.expr.data ) <- ma.data$genes$SystematicName

ma.flag.ok <- ma.data$other$gIsPosAndSignif == 1

raw.expr.data[ ! ma.flag.ok ] <- NA


# filter out control probes and probes with missing levels in any sample
# reorder samples
# change to log2 scale

ma.control.probe.bol <- ma.data$genes$ControlType != 0

ma.all.level.probe.bol <- rowSums( is.na( raw.expr.data ) ) == 0

real.expr.data <- log2( raw.expr.data[ 
    ! ma.control.probe.bol & ma.all.level.probe.bol, ma.sample ] )

if ( ! file.exists( data.dir ) )
    dir.create( data.dir )

save( real.expr.data, file = file.path( data.dir, "real_expr.data.rda" ) )

ma.probe <- rownames( real.expr.data )
ma.probe.num <- length( ma.probe )

rm( ma.data, ma.flag.ok, raw.expr.data )


# normalize real expression data and obtain offset

stdvec.condec.norm.real.result <- normalize.stdvec.condec( real.expr.data, 
    ma.sample.condition, verbose=TRUE )

stdvec.condec.norm.real.offset <- stdvec.condec.norm.real.result$offset

rm( stdvec.condec.norm.real.result )


# generate random data sets

random.seed <- 20140101
set.seed( random.seed )

expr.data.sample.mean <- rowMeans( real.expr.data )
expr.data.sample.sd <- apply( real.expr.data, 1, sd )

normal.expr.data <- 
    t( mapply( function( m, s ) rnorm( ma.sample.num, m, s ), 
        expr.data.sample.mean, expr.data.sample.sd ) )

lognormal.expr.data <- 
    t( mapply( function( m, s ) 
        m + s * scale( rlnorm( ma.sample.num, 0, 0.5 ) ), 
        expr.data.sample.mean, expr.data.sample.sd ) )

colnames( normal.expr.data ) <- ma.sample
colnames( lognormal.expr.data ) <- ma.sample

normal.expr.data <- sweep( normal.expr.data, 2, 
    stdvec.condec.norm.real.offset, "+" )
lognormal.expr.data <- sweep( lognormal.expr.data, 2, 
    stdvec.condec.norm.real.offset, "+" )

save( normal.expr.data, lognormal.expr.data, 
    file = file.path( data.dir, "random_expr.data.rda" ) )


# normalize real and random datasets with four methods

data.tag <- c( "real", "normal" )

norm.method <- c( "no", "median", "quantile", "median.condec", "stdvec.condec" )

if ( ! file.exists( boxplot.dir ) )
    dir.create( boxplot.dir )

for ( dtag in data.tag )
{
    expr.data <- get( paste0( dtag, ".expr.data" ) )
    
    p.value.graph.dir <- file.path( p.value.dir, dtag )
    
    # no normalization
    no.norm.data <- expr.data
    
    # median normalization
    median.offset <- apply( expr.data, 2, median )
    median.offset <- median.offset - mean( median.offset )
    median.norm.data <- sweep( expr.data, 2, median.offset )
    
    # quantile normalization
    quantile.norm.data <- normalizeBetweenArrays( expr.data, method="quantile" )
    
    # median condition-decomposition normalization
    median.condec.norm.result <- normalize.median.condec( expr.data, 
        ma.sample.condition, 0.5, verbose=TRUE )
    median.condec.norm.data <- median.condec.norm.result$data
    
    # standard-vector condition-decomposition normalization
    stdvec.condec.norm.result <- normalize.stdvec.condec( expr.data, 
        ma.sample.condition, verbose=TRUE, p.value.graph=p.value.graph.dir )
    stdvec.condec.norm.data <- stdvec.condec.norm.result$data
    
    # save normalization results    
    norm.data.var <- c( "no.norm.data", "median.norm.data", 
        "quantile.norm.data", "median.condec.norm.result", 
        "stdvec.condec.norm.result" )
    norm.file.name <- paste0( dtag, "_norm.data.rda" )
    save( list=norm.data.var, file = file.path( data.dir, norm.file.name ) )
    
    # plot expression levels
    y.at <- 4:10
    y.label <- as.character( y.at )
    y.label[ 2 * 1:3 ] <- ""
    
    for ( nmeth in norm.method )
    {
        norm.data <- get( paste0( nmeth, ".norm.data" ) )
        
        bp.file.name <- paste0( dtag, "_", nmeth, ".norm.pdf" )
        cairo_pdf( filename = file.path( boxplot.dir, bp.file.name ), 
            width=3.5, height=1.5 )
        
        par( mar = c( 0.1, 0.7, 0.1, 0.1 ), mgp = c( 0.75, 0, 0 ), tcl=-0.15 )
        
        boxplot( norm.data, xaxt="n", yaxt="n", 
            xlim = c( 4, ma.sample.num-3 ), ylim = c( 3.6, 10.4 ), 
            col=ma.sample.color, outline=FALSE, pars = list( boxwex=1, 
            medlwd=1.5, boxlty=0, whisklty=0 ) )
        
        axis( 2, at=y.at, labels=y.label , cex.axis=0.55 )
        
        dev.off()
    }  
}


# obtain graphs with convergence of standard vectors

stdvec.convergence.cond <- "Ag.NM300K.EC20.3d"

stdvec.convergence.condition <- ma.sample.condition
stdvec.convergence.condition[ ma.sample.condition != 
    stdvec.convergence.cond ] <- NA

data.tag <- c( "real", "normal", "lognormal" )

if ( ! file.exists( stdvec.dir ) )
    dir.create( stdvec.dir )

for ( dtag in data.tag )
{
    expr.data <- get( paste0( dtag, ".expr.data" ) )
    
    vector.graph.dir <- file.path( stdvec.dir, dtag )
    
    # plot convergence for one condition
    normalize.stdvec.condec( expr.data, stdvec.convergence.condition, 
        verbose=TRUE, vector.graph=vector.graph.dir )
}


# study between-condition variation, depending on probe selection

random.seed <- 20150201
set.seed( random.seed )

h0p.frac.num <- 300
h0p.frac.initial <- 1e-3
h0p.frac.final <- 1
h0.probe.frac <- h0p.frac.initial * 
    ( h0p.frac.final / h0p.frac.initial )^( 
        0 : ( h0p.frac.num - 1 ) / ( h0p.frac.num - 1 ) )
h0.probe.num <- floor( h0.probe.frac * ma.probe.num )
h0.probe.num <- h0.probe.num[ ! duplicated( h0.probe.num ) ]

data.tag <- c( "real", "normal" )

norm.method <- c( "median.condec", "stdvec.condec" )

between.condition.mean <- c( "grand.mean", "median.mean" )

if ( ! file.exists( bc.variation.dir ) )
    dir.create( bc.variation.dir )

for ( dtag in data.tag )
{
    norm.file.name <- paste0( dtag, "_norm.data.rda" )
    loaded <- load( file.path( data.dir, norm.file.name ) )
    
    stopifnot( c( "no.norm.data", "median.condec.norm.result", 
        "stdvec.condec.norm.result" ) %in% loaded )
    
    cond.grand.mean.sd <- data.frame( h0.probe.num )
    cond.median.mean.sd <- data.frame( h0.probe.num )
    
    for ( nmeth in norm.method )
    {        
        # obtain within-condition results
        norm.result <- get( paste0( nmeth, ".norm.result" ) )
        
        norm.within.cond.data <- sweep( no.norm.data, 2, 
            norm.result$within.condition.offset )
        
        # filter out zero signals (only needed in general, not for this dataset)
        norm.within.cond.data[ 
            norm.within.cond.data < max( norm.within.cond.data, na.rm=TRUE ) * 
                .Machine$double.eps ] <- NA
        
        # identify samples available per condition and per probe
        norm.within.cond.n <- sapply( ma.condition, function( cond ) 
            rowSums( ! is.na( norm.within.cond.data[ , 
                cond == ma.sample.condition ] ) ) )
        
        # calculate balanced within-condition means
        norm.bal.mean.n <- min( table( ma.sample.condition ) )
        
        # select probes that have at least bal.mean.n samples in all conditions
        norm.bal.probe.idx <- 
            which( rowSums( norm.within.cond.n < norm.bal.mean.n ) == 0 )
        
        norm.bal.within.cond.mean <- sapply( ma.condition, function( cond ) 
            if ( sum( cond == ma.sample.condition ) == norm.bal.mean.n )
                rowMeans( norm.within.cond.data[ norm.bal.probe.idx, 
                    cond == ma.sample.condition ] )
            else
                rowMeans( t( apply( 
                    t( norm.within.cond.data[ norm.bal.probe.idx, 
                        cond == ma.sample.condition ] ), 
                    2, function( ed ) sample( na.omit( ed ), norm.bal.mean.n ) 
                ) ) )
        )
        
        # for F-statistic
        # select probes that have at least 2 samples in all conditions
        norm.probe.idx <- 
            which( rowSums( norm.within.cond.n < 2 ) == 0 )
        
        norm.within.cond.n <- norm.within.cond.n[ norm.probe.idx, ]
        
        # calculate within-condition means
        norm.within.cond.mean <- sapply( ma.condition, function( cond ) 
            rowMeans( norm.within.cond.data[ norm.probe.idx, 
                cond == ma.sample.condition ], 
                na.rm=TRUE ) )
        norm.within.cond.grand.mean <- colMeans( norm.within.cond.mean )
        
        # calculate within-condition variances
        norm.within.cond.var <- sapply( ma.condition, function( cond ) 
            apply( norm.within.cond.data[ norm.probe.idx, 
                cond == ma.sample.condition ], 1, var, na.rm=TRUE )
        ) * ( norm.within.cond.n - 1 )
        norm.within.cond.var <- rowSums( norm.within.cond.var ) / 
                ( rowSums( norm.within.cond.n ) - ma.condition.num )
        
        # calculate within-condition medians
        norm.within.cond.median <- apply( norm.within.cond.data, 2, 
            median, na.rm=TRUE )
        norm.within.cond.median.mean <- sapply( ma.condition, function( cond ) 
            mean( norm.within.cond.median[ cond == ma.sample.condition ] ) )
        
        for ( h0p.random in c( TRUE, FALSE ) )
        {
            norm.h0p.method <- paste0( nmeth, 
                ifelse( h0p.random, ".random", ".pvalue" ) )
            
            for ( h0p.num in h0.probe.num )
            {
                cat( paste( dtag, norm.h0p.method, h0p.num, "\n" ) )
                
                norm.between.cond.result <- switch( nmeth, 
                    stdvec.condec = normalize.standard.vector.probe( 
                        norm.bal.within.cond.mean, norm.within.cond.mean, 
                        norm.within.cond.var, norm.within.cond.n, 
                        "between.condition", h0p.random, h0p.num, verbose=TRUE, 
                        vector.graph=NULL ), 
                    median.condec = normalize.median.selection.probe( 
                        norm.bal.within.cond.mean, norm.within.cond.mean, 
                        norm.within.cond.var, norm.within.cond.n, 0.5, 
                        h0p.random, h0p.num, verbose=TRUE ) )
                
                norm.between.cond.grand.mean.sd <- 
                    sd( norm.within.cond.grand.mean - 
                        norm.between.cond.result$offset )
                
                norm.between.cond.median.mean.sd <- 
                    sd( norm.within.cond.median.mean - 
                        norm.between.cond.result$offset )
                
                cond.grand.mean.sd[[ norm.h0p.method ]][ 
                    cond.grand.mean.sd$h0.probe.num == h0p.num ] <- 
                    norm.between.cond.grand.mean.sd
                
                cond.median.mean.sd[[ norm.h0p.method ]][ 
                    cond.median.mean.sd$h0.probe.num == h0p.num ] <- 
                    norm.between.cond.median.mean.sd
            }
        }
    }
    
    # save variance of grand means
    grand.mean.sd.file.name <- paste0( dtag, "_grand.mean.sd.rda" )
    save( cond.grand.mean.sd, 
        file = file.path( data.dir, grand.mean.sd.file.name ) )

    # save variance of median means
    median.mean.sd.file.name <- paste0( dtag, "_median.mean.sd.rda" )
    save( cond.median.mean.sd, 
        file = file.path( data.dir, median.mean.sd.file.name ) )
    
    # plot variance of grand means and median means
    for ( bcmean in between.condition.mean )
    {
        bcmp.file.name <- paste0( dtag, "_", bcmean, ".sd.pdf" )
        cairo_pdf( filename = file.path( bc.variation.dir, bcmp.file.name ), 
            width=3.5, height=2.8 )
        
        par( mar = c( 1.7, 1.45, 0.15, 0.15 ), mgp = c( 0.75, 0, 0 ), 
            tcl=-0.15 )
        
        bcmp.xlim <- range( h0.probe.num )
        x.tick.pos <- c( 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000 )
        x.tick.lab <- c( "20", "", "100", "", "", "1000", "", "", "10000", "" )
        
        if ( bcmean == "grand.mean" )
        {
            if ( dtag == "real" ) {
                bcmp.ylim <- 0.5 * 10^c( -1.7, 0.1 )
                y.tick.pos <- c( 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 )
                y.tick.lab <- c( "0.01", "", "", "0.1", "", "0.5" )
            }
            else {    # dtag == "normal"
                bcmp.ylim <- 0.5 * 10^c( -2.7, 0.1 )
                y.tick.pos <- 
                    c( 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 )
                y.tick.lab <- 
                    c( "0.001", "", "", "0.01", "", "", "0.1", "", "0.5" )
            }
        }
        else    # bcmean == "median.mean"
        {
            if ( dtag == "real" )
                bcmp.ylim <- 0.5 * 10^c( -1.8, 0.1 )
            else    # dtag == "normal"
                bcmp.ylim <- 0.5 * 10^c( -1.9, 0.1 )
            y.tick.pos <- c( 0.01, 0.02, 0.05, 0.1, 0.2, 0.5 )
            y.tick.lab <- c( "0.01", "", "", "0.1", "", "0.5" )
        }
        
        norm.method.col <- c( "black", "blue" )
        names( norm.method.col ) <- norm.method
        
        plot( 1, type="n", cex.lab=0.55, log="xy", xaxt="n", yaxt="n", 
            xlim=bcmp.xlim , ylim=bcmp.ylim, 
            xlab="# Gene probes in between-condition normalization", 
            ylab="Between-condition variation" )
        
        axis( 1, at=x.tick.pos, labels=x.tick.lab, cex.axis=0.55 )
        axis( 2, at=y.tick.pos, labels=y.tick.lab, cex.axis=0.55 )
        
        segments( 20, 0.5, 20e3, 0.5*10^-1.5, lty=2, lwd=2 )
        
        for ( nmeth in norm.method )
        {
            norm.result <- get( paste0( nmeth, ".norm.result" ) )
            
            norm.h0.probe.num <- 
                length( norm.result$between.condition.h0.probe )
            
            if ( bcmean == "grand.mean" ) {
                norm.cond.mean <- sapply( ma.condition, function( cond ) 
                    rowMeans( norm.result$data[ , 
                        cond == ma.sample.condition ], na.rm=TRUE ) )
                norm.cond.bcmean.sd <- sd( colMeans( norm.cond.mean ) )
            }
            else {    # bcmean == "median.mean"
                norm.cond.median <- apply( norm.result$data, 2, median, 
                    na.rm=TRUE )
                norm.cond.median.mean <- sapply( ma.condition, function( cond ) 
                    mean( norm.cond.median[ cond == ma.sample.condition ] ) )
                norm.cond.bcmean.sd <- sd( norm.cond.median.mean )
            }
            
            points( norm.h0.probe.num, norm.cond.bcmean.sd, 
                col = norm.method.col[ nmeth ], pch=1, cex=2, lwd=2 )
            
            cond.bcmean.sd <- get( paste0( "cond.", bcmean, ".sd" ) )
            
            points( cond.bcmean.sd$h0.probe.num, 
                cond.bcmean.sd[[ paste0( nmeth, ".random" ) ]], cex=0.8, 
                col = norm.method.col[ nmeth ], pch=1 )
            
            points( cond.bcmean.sd$h0.probe.num, 
                cond.bcmean.sd[[ paste0( nmeth, ".pvalue" ) ]], cex=0.6, 
                col = norm.method.col[ nmeth ], pch=20 )
        }
        
        dev.off()
    }
}


# study the number of differentially expressed genes in treatment-vs-control 
# comparisons

ma.treatment.test <- as.vector( sapply( ma.treatment, function( mt )
    paste0( mt, c( ".fc", ".pv", ".fdr" ) ) ) )

dtag <- "real"

norm.method <- c( "median", "quantile", "median.condec", "stdvec.condec" )

norm.data.var <- c( "median.norm.data", "quantile.norm.data", 
    "median.condec.norm.result", "stdvec.condec.norm.result" )

norm.file.name <- paste0( dtag, "_norm.data.rda" )
loaded <- load( file.path( data.dir, norm.file.name ) )
stopifnot( norm.data.var %in% loaded )

median.condec.norm.data <- median.condec.norm.result$data
stdvec.condec.norm.data <- stdvec.condec.norm.result$data

diff.expr.method <- c( "ttest", "limma" )

alpha <- 0.05

if ( ! file.exists( diff.expr.dir ) )
    dir.create( diff.expr.dir )

# build limma model matrix

condition.factor <- factor( ma.sample.condition, levels=ma.condition )

model.design <- model.matrix( ~0 + condition.factor )
rownames( model.design ) <- ma.sample
colnames( model.design ) <- ma.condition

diff.expr.contrast <- mapply( function( t, c ) {
    de.cont <- rep( 0, ma.condition.num )
    de.cont[ ma.condition == t ] <- 1
    de.cont[ ma.condition == c ] <- -1
    de.cont
}, ma.treatment, ma.control )
rownames( diff.expr.contrast ) <- ma.condition

# iterate over stat methods and then over normalizations

for ( demeth in diff.expr.method )
{
    for ( nmeth in norm.method )
    {
        norm.data <- get( paste0( nmeth, ".norm.data" ) )
        
        stopifnot( rownames( norm.data ) == ma.probe )
        stopifnot( colnames( norm.data ) == ma.sample )
        
        if ( demeth == "ttest" )
        {
            diff.expr <-  mapply( function( t, c ) {
                tc.col <- ma.sample.condition == t | ma.sample.condition == c
                tc.factor <- factor( ( 1*( ma.sample.condition == t ) + 
                        2*( ma.sample.condition == c ) )[ tc.col ] )
                rtt <- rowttests( norm.data[ , tc.col ], tc.factor )
                rtt$fdr <- p.adjust( rtt$p.value, method="fdr" )
                c( rtt$dm, rtt$p.value, rtt$fdr )
            }, ma.treatment, ma.control )
        }
        else    # demeth == "limma"
        {
            model.fit <- lmFit( norm.data, model.design )
            diff.expr.fit <- contrasts.fit( model.fit, diff.expr.contrast )
            diff.expr.efit <- eBayes( diff.expr.fit )
            diff.expr <- sapply( ma.treatment, function( mt ) {
                de.table <- topTable( diff.expr.efit, coef=mt, number=Inf, 
                    sort.by="none" )
                c( de.table$logFC, de.table$P.Value, de.table$adj.P.Val )
            } )
        }
        
        dim( diff.expr ) <- c( ma.probe.num, length( ma.treatment.test ) )
        rownames( diff.expr ) <- ma.probe
        colnames( diff.expr ) <- ma.treatment.test
        
        assign( paste0( nmeth, ".diff.expr" ), 
            data.frame( probe=ma.probe, diff.expr, stringsAsFactors=FALSE ) )
    }
    
    diff.expr.file.name <- paste0( dtag, "_", demeth, ".diff.expr.rda" )
    save( list = paste0( norm.method, ".diff.expr" ), 
        file = file.path( data.dir, diff.expr.file.name ) )
    
    # find the number of differentially expressed genes
    
    for ( nmeth in norm.method )
    {
        diff.expr <- get( paste0( nmeth, ".diff.expr" ) )
        diff.expr.fdr <- diff.expr[ , grep( "\\.fdr$", names( diff.expr ) ) ]
        diff.expr.num <- colSums( diff.expr.fdr < alpha, na.rm=TRUE )
        
        names( diff.expr.num ) <- sub( "\\.fdr$", "", names( diff.expr.num ) )
        stopifnot( names( diff.expr.num ) == ma.treatment )
        
        assign( paste0( nmeth, ".diff.expr.num" ), diff.expr.num )
    }
    
    diff.expr.num.file.name <- paste0( dtag, "_", demeth, ".diff.expr.num.rda" )
    save( list = paste0( norm.method, ".diff.expr.num" ), 
        file = file.path( data.dir, diff.expr.num.file.name ) )
    
    # plot the number of differentially expressed genes
    
    dep.file.name <- paste0( dtag, "_", demeth, ".diff.expr.pdf" )
    cairo_pdf( filename = file.path( diff.expr.dir, dep.file.name ), width=3.5, 
        height=2.8 )
    
    par( mar = c( 1.7, 1.45, 0.15, 0.15 ), mgp = c( 0.75, 0, 0 ), tcl=-0.15 )
    
    norm.method.pch <- c( 1, 2, 19, 17 )
    norm.method.col <- c( "black", "red3", "green3", "blue3" )
    
    names( norm.method.pch ) <- norm.method
    names( norm.method.col ) <- norm.method
    
    plot( 1, type="n", cex.lab=0.55, cex.axis=0.55, log="y", xaxt="n", 
        xlim = c( 1, ma.treatment.num ), ylim = c( 1, ma.probe.num ), 
        xlab="Treatment vs control comparison", 
        ylab="# Differentially expressed gene probes" )
    
    axis( 1, at = seq( 5, ma.treatment.num, 5 ), cex.axis=0.55 )
    
    diff.expr.num.order <- order( stdvec.condec.diff.expr.num )
    
    for ( nmeth in norm.method )
    {
        diff.expr.num <- get( paste0( nmeth, ".diff.expr.num" ) )
        diff.expr.num.plot <- diff.expr.num[ diff.expr.num.order ]
        
        points( ( 1 : ma.treatment.num )[ diff.expr.num.plot > 0 ], 
            diff.expr.num.plot[ diff.expr.num.plot > 0 ], cex=1, lwd=1, 
            pch=norm.method.pch[ nmeth ], col=norm.method.col[ nmeth ] )
    }
    
    dev.off()
}


# study the magnitude of differential expression

dtag <- "real"

demeth <- "limma"

norm.method <- c( "quantile", "stdvec.condec" )

diff.expr.var <- paste0( norm.method, ".diff.expr" )
diff.expr.file.name <- paste0( dtag, "_", demeth, ".diff.expr.rda" )
loaded <- load( file.path( data.dir, diff.expr.file.name ) )
stopifnot( diff.expr.var %in% loaded )

diff.expr.num.var <- paste0( norm.method, ".diff.expr.num" )
diff.expr.num.file.name <- paste0( dtag, "_", demeth, ".diff.expr.num.rda" )
loaded <- load( file.path( data.dir, diff.expr.num.file.name ) )
stopifnot( diff.expr.num.var %in% loaded )

diff.expr.num.order <- order( stdvec.condec.diff.expr.num )

# write table of sorted treatment-vs-control comparisons
treat.cont.file.name <- paste0( dtag, "_", demeth, ".treatment.control.csv" )
write.csv( cbind( 1:ma.treatment.num, 
    ma.treatment[ diff.expr.num.order ], 
    ma.control[ diff.expr.num.order ], 
    stdvec.condec.diff.expr.num[ diff.expr.num.order ] ), 
    file.path( table.dir, treat.cont.file.name ), row.names=FALSE )

for ( nmeth in norm.method )
{
    diff.expr.fc <- lapply( ma.treatment, function( mt ) {
        de.data <- get( paste0( nmeth, ".diff.expr" ) )
        de.fdr <- de.data[[ paste0( mt, ".fdr" ) ]]
        de.fc <- de.data[[ paste0( mt, ".fc" ) ]][ de.fdr < alpha ]
        if ( is.null( de.fc ) ) de.fc else abs( de.fc )
    } )
    
    # plot fold change of differentially expressed genes
    
    fcp.file.name <- paste0( dtag, "_", demeth, ".diff.expr.fc_", nmeth, 
        ".pdf" )
    cairo_pdf( filename = file.path( diff.expr.dir, fcp.file.name ), width=3.5, 
        height=2.8 )
    
    par( mar = c( 1.7, 1.45, 0.15, 0.15 ), mgp = c( 0.75, 0, 0 ), tcl=-0.15 )
    
    boxplot( diff.expr.fc[ diff.expr.num.order ], cex.lab=0.55, cex.axis=0.55, 
        log="y", xaxt="n", 
        xlim = c( 1, ma.treatment.num ), ylim = c( 0.3, 11 ), 
        xlab="Treatment vs control comparison", ylab="| DEGP fold change |", 
        col = ma.treatment.color[ diff.expr.num.order ], pch=20, cex=0.5, 
        lwd=0.5 )
    
    axis( 1, at = seq( 5, ma.treatment.num, 5 ), cex.axis=0.55 )
    
    abline( h=log2( 1.5 ), lty=2, lwd=2 )
    abline( h=1, lty=2, lwd=2 )
    
    dev.off()
}

