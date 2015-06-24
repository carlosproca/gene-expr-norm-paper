# Copyright (c) 2015, Universitat Rovira i Virgili (Spain), University of 
# Aveiro (Portugal) & Aarhus University (Denmark)
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


# Adapts median condition-decompositon normalization for study across probes


normalize.median.selection.probe <- function( edata, edata.fstat, 
    within.cond.var, within.cond.n, norm.prob, h0.probe.random, h0.probe.num, 
    verbose )
{
    iter.max <- 100
    median.offset.accum.step.max <- 10
    median.offset.accum.threshold <- 0.8    # 0.1 -> 0.8 increase needed for 
    median.offset.single.threshold <- 0.01  # convergence with very small 
                                            # h0.probe.num
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    last.median.offset <- vector( "numeric" )
    last.median.h0.probe <- vector( "character" )
    
    norm.median.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.median.offset.ratio <- vector( "numeric" )
    norm.median.offset.accum.step <- vector( "numeric" )
    norm.median.h0.probe.num <- vector( "numeric" )

    median.offset <- rep( 0, ncol( edata ) )
    
    if ( h0.probe.random )
        h0.probe <- sample( rownames( edata.fstat ), h0.probe.num )
    else
        h0.probe <- NULL
    
    iter <- 0
    median.offset.accum.step <- 0
    median.offset.ratio <- 1
    
    while ( iter < iter.max && 
        median.offset.accum.step < median.offset.accum.step.max &&
        median.offset.ratio >= median.offset.single.threshold )
    {
        iter <- iter + 1
        
        # obtain next step of median offset     
        median.offset.step <- calculate.median.offset.probe( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, norm.prob, 
            h0.probe, h0.probe.num, iter )
        
        median.offset.delta <- median.offset.step$value
        median.h0.probe <- median.offset.step$h0.probe
        
        median.h0.probe.num <- length( median.h0.probe )
        
        # check errors
        if ( any( is.nan( median.offset.delta ) ) )
            stop( "NaN error in normalize.median.condec.probe" )
        
        # update total median offset
        median.offset <- median.offset + median.offset.delta
        median.offset <- median.offset - mean( median.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, median.offset )
        if ( ! is.null( edata.fstat ) )
            last.norm.data.fstat <- sweep( edata.fstat, 2, median.offset )
        
        # check convergence
        median.offset.sd <- sd( median.offset )
        median.offset.delta.sd <- sd( median.offset.delta )
        
        median.offset.ratio <- median.offset.delta.sd / median.offset.sd
        
        if ( median.offset.ratio < median.offset.accum.threshold )
            median.offset.accum.step <- median.offset.accum.step + 1
        else
            median.offset.accum.step <- 0
        
        # store last results
        last.median.offset <- median.offset
        last.median.h0.probe <- median.h0.probe
        
        # store step results
        norm.median.offset <- rbind( norm.median.offset, median.offset )
        norm.median.offset.ratio <- c( norm.median.offset.ratio, 
            median.offset.ratio )
        norm.median.offset.accum.step <- c( norm.median.offset.accum.step, 
            median.offset.accum.step )
        norm.median.h0.probe.num <- c( norm.median.h0.probe.num, 
            median.h0.probe.num )
        
        if ( verbose )
            cat( sprintf( "  %2d %g %g %02d %d\n", iter, median.offset.sd, 
                median.offset.ratio, median.offset.accum.step, 
                median.h0.probe.num ) )
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( iter == iter.max )
        stop( "no convergence in normalize.median.condec.probe" )
    
    dimnames( norm.median.offset ) <- NULL
    
    norm.median.convergence <- list( offset = norm.median.offset, 
        offset.ratio = norm.median.offset.ratio, 
        offset.accum.step = norm.median.offset.accum.step, 
        h0.probe.num = norm.median.h0.probe.num )
    
    list( offset = last.median.offset, h0.probe = last.median.h0.probe, 
        convergence = norm.median.convergence )
}


calculate.median.offset.probe <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, h0.probe, h0.probe.num, iter )
{
    if ( is.null( h0.probe ) )
    {
        # calculate f statistics for each probe
        expr.k <- ncol( within.cond.n )
        expr.n <- rowSums( within.cond.n )
        
        expr.grand.mean <- rowSums( edata.fstat * within.cond.n ) / expr.n
        
        between.cond.var <- rowSums( ( edata.fstat - expr.grand.mean )^2 * 
                within.cond.n ) / ( expr.k - 1 )
        
        expr.f <- between.cond.var / within.cond.var
        
        expr.p.value <- pf( expr.f, df1 = expr.k-1, df2 = expr.n-expr.k, 
            lower.tail = FALSE )
        
        h0.probe.idx <- order( expr.p.value, decreasing=TRUE )[ 
            1 : h0.probe.num ]
        h0.probe <- names( expr.p.value )[ h0.probe.idx ]
    }
    
    # identify probes for normalization
    median.probe <- intersect( h0.probe, rownames( edata ) )
    
    # find median offset using only H0 probes
    median.offset <- apply( edata[ median.probe, ], 2, quantile, 
        probs=norm.prob )
    median.offset <- median.offset - mean( median.offset )
    
    list( value = median.offset, h0.probe = h0.probe )
}

