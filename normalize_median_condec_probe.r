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


# Adapts median condition-decompositon normalization for study across gene 
# probes


normalize.median.selection.probe <- function( edata, edata.fstat, 
    within.cond.var, within.cond.n, norm.prob, h0.probe.random, h0.probe.num, 
    verbose )
{
    if ( h0.probe.random )
    {
        median.h0.probe <- sample( rownames( edata ), h0.probe.num )
        
        median.offset <- apply( edata[ median.h0.probe, ], 2, quantile, 
            probs=norm.prob )
        median.offset <- median.offset - mean( median.offset )
        
        if ( verbose )
            cat( "\n" )
        
        norm.median.result <- list( offset = median.offset, 
            h0.probe = median.h0.probe )
        
        return( norm.median.result )
    }
    
    edata.probe.num <- nrow( edata )
    
    iter.max <- 100
    single.threshold <- 0.001
    accum.threshold <- ifelse( h0.probe.num < 0.05*edata.probe.num, 1.0, 0.1 )
    offset.accum.step.threshold <- 10
    
    median.offset <- rep( 0, ncol( edata ) )
    
    norm.median.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.median.offset.sd <- vector( "numeric" )
    norm.median.offset.delta.sd <- vector( "numeric" )
    norm.median.offset.accum.step <- vector( "numeric" )
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    median.offset.accum.step <- 0
    
    iter <- 0
    overall.convergence <- FALSE
    
    while ( iter < iter.max && ! overall.convergence )
    {
        iter <- iter + 1
        
        # obtain next step of median offset     
        median.offset.step <- calculate.median.offset.probe( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, norm.prob, 
            h0.probe.num, iter )
        
        median.offset.delta <- median.offset.step$value
        
        median.h0.probe <- median.offset.step$h0.probe
        
        # check errors
        if ( any( is.nan( median.offset.delta ) ) )
            stop( "NaN error in normalize.median.condec.probe" )
        
        # update total median offset
        median.offset <- median.offset + median.offset.delta
        median.offset <- median.offset - mean( median.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, median.offset )
        last.norm.data.fstat <- sweep( edata.fstat, 2, median.offset )
        
        # check convergence
        median.offset.sd <- sd( median.offset )
        median.offset.delta.sd <- sd( median.offset.delta )
        
        median.offset.delta.sd.ratio <- ifelse( median.offset.sd > 0, 
            median.offset.delta.sd / median.offset.sd, 1 )
        
        if ( median.offset.delta.sd == 0 )
            median.offset.ratio <- 0
        else if ( median.offset.sd == 0 )
            median.offset.ratio <- 1
        else
            median.offset.ratio <- median.offset.delta.sd / median.offset.sd
        
        median.offset.accum.step <- ifelse( 
            median.offset.ratio < accum.threshold, 
            median.offset.accum.step + 1, 0 )
        
        median.convergence <- median.offset.ratio < single.threshold || 
            median.offset.accum.step > offset.accum.step.threshold
        
        overall.convergence <- median.convergence
        
        # store last results
        last.median.offset <- median.offset
        
        # store step results
        norm.median.offset <- rbind( norm.median.offset, median.offset )
        norm.median.offset.sd <- c( norm.median.offset.sd, median.offset.sd )
        norm.median.offset.delta.sd <- c( norm.median.offset.delta.sd, 
            median.offset.delta.sd )
        norm.median.offset.accum.step <- c( norm.median.offset.accum.step, 
            median.offset.accum.step )
        
        if ( verbose )
            cat( sprintf( "  %2d %g %g %02d\n", iter, 
                median.offset.sd, median.offset.delta.sd.ratio, 
                median.offset.accum.step ) )
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( ! overall.convergence )
        stop( "No convergence in normalize.median.condec.probe" )
    
    # remove sample or condition names from step results
    dimnames( norm.median.offset ) <- NULL
    
    norm.median.convergence <- list( offset = norm.median.offset, 
        offset.sd = norm.median.offset.sd, 
        offset.delta.sd = norm.median.offset.delta.sd, 
        offset.accum.step = norm.median.offset.accum.step )
    
    norm.median.result <- list( offset = last.median.offset, 
        convergence = norm.median.convergence, 
        h0.probe = median.h0.probe )
    
    norm.median.result
}


calculate.median.offset.probe <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, h0.probe.num, iter )
{
    # calculate f statistics for each probe
    expr.k <- length( within.cond.n )
    expr.n <- sum( within.cond.n )
    
    expr.grand.mean <- apply( edata.fstat, 1, function( ef ) 
        sum( ef * within.cond.n ) ) / expr.n
    
    between.cond.var <- apply( ( edata.fstat - expr.grand.mean )^2, 1, 
        function( ef2 ) sum( ef2 * within.cond.n ) ) / ( expr.k - 1 )
    
    expr.f <- between.cond.var / within.cond.var
    
    expr.f <- na.omit( expr.f )  # in case of 0/0
    attr( expr.f, "na.action" ) <- NULL
    
    expr.p.value <- pf( expr.f, df1 = expr.k-1, df2 = expr.n-expr.k, 
        lower.tail = FALSE )
    
    h0.probe.idx <- order( expr.p.value, decreasing=TRUE )[ 
        1 : h0.probe.num ]
    h0.probe <- names( expr.p.value )[ h0.probe.idx ]
    
    # calculate offset
    median.offset <- apply( edata[ h0.probe, ], 2, quantile, probs=norm.prob )
    median.offset <- median.offset - mean( median.offset )
    
    median.result <- list( value = median.offset, h0.probe = h0.probe )
    
    median.result
}

