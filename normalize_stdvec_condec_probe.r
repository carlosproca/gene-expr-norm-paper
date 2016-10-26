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


# Adapts standard-vector condition-decompositon normalization for study across 
# gene probes


normalize.standard.vector.probe <- function( edata, edata.fstat, 
    within.cond.var, within.cond.n, h0.probe.random, h0.probe.num, verbose )
{
    edata.probe.num <- nrow( edata )
    
    iter.max <- 100
    single.threshold <- 0.01
    accum.threshold <- ifelse( h0.probe.num < 0.05*edata.probe.num, 1.0, 0.1 )
    offset.accum.step.threshold <- 10
    
    stdvec.offset <- rep( 0, ncol( edata ) )
    
    norm.stdvec.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.stdvec.offset.sd <- vector( "numeric" )
    norm.stdvec.offset.stderr <- vector( "numeric" )
    norm.stdvec.offset.delta.sd <- vector( "numeric" )
    norm.stdvec.offset.accum.step <- vector( "numeric" )
    norm.stdvec.numerical.demand <- vector( "numeric" )
    
    if ( ncol( edata ) >= 3 )
        norm.stdvec.watson.u2 <- matrix( nrow=0, 
            ncol = ( ncol( edata ) - 1 ) %/% 3 + 1 )
    else
        norm.stdvec.watson.u2 <- NULL
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    stdvec.offset.accum.step <- 0
    
    if ( h0.probe.random )
        h0.probe <- sample( rownames( edata ), h0.probe.num )
    else
        h0.probe <- NULL
    
    iter <- 0
    overall.convergence <- FALSE
    
    while ( iter < iter.max && ! overall.convergence )
    {
        iter <- iter + 1
        
        # obtain next step of standard vector offset
        stdvec.offset.step <- calculate.stdvec.offset.probe( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, h0.probe, 
            h0.probe.num, iter )
        
        stdvec.offset.delta <- stdvec.offset.step$value
        stdvec.offset.stderr <- stdvec.offset.step$stderr
        stdvec.numerical.demand <- stdvec.offset.step$numerical.demand
        stdvec.watson.u2 <- stdvec.offset.step$watson.u2
        
        stdvec.h0.probe <- stdvec.offset.step$h0.probe
        
        # check errors
        if ( any( is.nan( stdvec.offset.delta ) ) || 
                is.nan( stdvec.offset.stderr ) )
            stop( "NaN error in normalize.stdvec.condec.probe" )
        
        if ( stdvec.numerical.demand < .Machine$double.eps * 10^3 )
            stop( "Numerical error in normalize.stdvec.condec.probe" )
        
        # update total standard vector offset
        stdvec.offset <- stdvec.offset + stdvec.offset.delta
        stdvec.offset <- stdvec.offset - mean( stdvec.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, stdvec.offset )
        last.norm.data.fstat <- sweep( edata.fstat, 2, stdvec.offset )
        
        # check convergence
        stdvec.offset.sd <- sd( stdvec.offset )
        stdvec.offset.delta.sd <- sd( stdvec.offset.delta )
        
        stdvec.offset.stderr.ratio <- ifelse( stdvec.offset.sd > 0, 
            stdvec.offset.stderr / stdvec.offset.sd, 1 )
        stdvec.offset.delta.sd.ratio <- ifelse( stdvec.offset.sd > 0, 
            stdvec.offset.delta.sd / stdvec.offset.sd, 1 )
        
        if ( stdvec.offset.delta.sd == 0 )
            stdvec.offset.ratio <- 0
        else if ( stdvec.offset.stderr == 0 )
            stdvec.offset.ratio <- 1
        else
            stdvec.offset.ratio <- stdvec.offset.delta.sd / stdvec.offset.stderr
        
        stdvec.offset.accum.step <- ifelse( 
            stdvec.offset.ratio < accum.threshold, 
            stdvec.offset.accum.step + 1, 0 )
        
        stdvec.convergence <- stdvec.offset.ratio < single.threshold || 
            stdvec.offset.accum.step > offset.accum.step.threshold
        
        overall.convergence <- stdvec.convergence
        
        # store last results
        last.stdvec.offset <- stdvec.offset
        
        # store step results
        norm.stdvec.offset <- rbind( norm.stdvec.offset, stdvec.offset )
        norm.stdvec.offset.sd <- c( norm.stdvec.offset.sd, stdvec.offset.sd )
        norm.stdvec.offset.stderr <- c( norm.stdvec.offset.stderr, 
            stdvec.offset.stderr )
        norm.stdvec.offset.delta.sd <- c( norm.stdvec.offset.delta.sd, 
            stdvec.offset.delta.sd )
        norm.stdvec.offset.accum.step <- c( norm.stdvec.offset.accum.step, 
            stdvec.offset.accum.step )
        norm.stdvec.numerical.demand <- c( norm.stdvec.numerical.demand, 
            stdvec.numerical.demand )
        norm.stdvec.watson.u2 <- rbind( norm.stdvec.watson.u2, 
            stdvec.watson.u2 )
        
        if ( verbose )
        {
            if ( ! is.null( stdvec.watson.u2 ) )
                stdvec.watson.u2.char <- paste0( signif( stdvec.watson.u2, 6 ), 
                    collapse=" " )
            else
                stdvec.watson.u2.char <- ""
            
            cat( sprintf( "  %2d %g %g %g %02d %g [%s]\n", iter, 
                stdvec.offset.sd, stdvec.offset.stderr.ratio, 
                stdvec.offset.delta.sd.ratio, stdvec.offset.accum.step, 
                stdvec.numerical.demand, stdvec.watson.u2.char ) )
        }
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( ! overall.convergence )
        stop( "No convergence in normalize.stdvec.condec.probe" )
    
    # remove sample or condition names from step results
    dimnames( norm.stdvec.offset ) <- NULL
    if ( ! is.null( norm.stdvec.watson.u2 ) )
        dimnames( norm.stdvec.watson.u2 ) <- NULL
    
    norm.stdvec.convergence <- list( offset = norm.stdvec.offset, 
        offset.sd = norm.stdvec.offset.sd, 
        offset.stderr = norm.stdvec.offset.stderr, 
        offset.delta.sd = norm.stdvec.offset.delta.sd, 
        offset.accum.step = norm.stdvec.offset.accum.step, 
        numerical.demand = norm.stdvec.numerical.demand, 
        watson.u2 = norm.stdvec.watson.u2 )
    
    norm.stdvec.result <- list( offset = last.stdvec.offset, 
        convergence = norm.stdvec.convergence, 
        h0.probe = stdvec.h0.probe )
    
    norm.stdvec.result
}


calculate.stdvec.offset.probe <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, h0.probe, h0.probe.num, iter )
{
    stdvec.trim <- 0.01
    
    if ( is.null( h0.probe ) )
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
    }
    
    # identify probes for normalization
    expr.var <- apply( edata.fstat[ h0.probe, ], 1, var )
    expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
    
    stdvec.probe <- h0.probe[ 
        expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
        expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ]
    
    # center and scale expression data
    expr.mean <- rowMeans( edata[ stdvec.probe, ] )
    expr.centered <- sweep( edata[ stdvec.probe, ], 1, expr.mean, "-" )
    
    expr.sd.inv <- 1 / apply( expr.centered, 1, sd )
    expr.scaled <- sweep( expr.centered, 1, expr.sd.inv, "*" )
    
    # calculate offset
    expr.sd.inv.sum <- sum( expr.sd.inv )
    stdvec.offset <- apply( expr.scaled, 2, sum ) / expr.sd.inv.sum
    stdvec.offset <- stdvec.offset - mean( stdvec.offset )
    
    # estimate error and numerical demand
    stdvec.offset.stderr <- sqrt( length( stdvec.probe ) ) / expr.sd.inv.sum
    expr.sd.inv.min <- min( expr.sd.inv )
    stdvec.numerical.demand <- expr.sd.inv.min / expr.sd.inv.sum
    
    # calculate density distribution of standard vector angles
    
    dimension.num <- ncol( expr.scaled )
    
    if ( dimension.num < 3 )
        theta.watson.u2 <- NULL
    else
    {
        # identify condition groups
        dimension.group.num <- ( dimension.num - 1 ) %/% 3 + 1
        
        dimension.group <- lapply( 1 : dimension.group.num, function ( g )
            if ( g < dimension.group.num )
                ( 3*g - 2 ) : ( 3*g )
            else
                ( dimension.num - 2 ) : dimension.num )
        
        theta.watson.u2 <- vector( "numeric", dimension.group.num )
        
        uv <- matrix( c( 0, -1/sqrt(2), 1/sqrt(2), 
            2/sqrt(6), -1/sqrt(6), -1/sqrt(6) ), nrow=3 )
        
        for ( dim.group.idx in 1:dimension.group.num )
        {
            # select expression values for each condition group
            expr.dim.group <- expr.scaled[ , 
                dimension.group[[ dim.group.idx ]] ]
            
            if ( dimension.group.num > 1 )
            {
                # re-standardize again for this group
                expr.dim.group.mean <- rowMeans( expr.dim.group )
                expr.dim.group.centered <- sweep( expr.dim.group, 1, 
                    expr.dim.group.mean, "-" )
                
                expr.dim.group.sd <- apply( expr.dim.group.centered, 1, sd )
                expr.dim.group.sel <- expr.dim.group.sd != 0
                expr.dim.group <- sweep( 
                    expr.dim.group.centered[ expr.dim.group.sel, ], 1, 
                    expr.dim.group.sd[ expr.dim.group.sel ], "/" )
            }
            
            expr.uv <- expr.dim.group %*% uv
            expr.u <- expr.uv[ , 1 ]
            expr.v <- expr.uv[ , 2 ]
            
            # calculate density distribution of angles
            theta.density.n <- 2^11
            theta.density.adjust <- 0.5
            
            expr.theta <- atan2( expr.v, expr.u )
            expr.theta.ex <- c( expr.theta, 
                expr.theta + ifelse( expr.theta > 0, -2*pi, 2*pi ) )
            
            expr.theta.density <- density( expr.theta.ex, 
                adjust=theta.density.adjust, n=theta.density.n )
            expr.theta.density.sel <- expr.theta.density$x > -pi & 
                expr.theta.density$x <= pi
            
            # calculate density distribution of angles after permutations
            expr.theta.permu <- cbind( expr.theta.ex, 
                expr.theta.ex + (2*pi)/3, 
                expr.theta.ex - (2*pi)/3, 
                - expr.theta.ex + pi, 
                - expr.theta.ex + pi/3, 
                - expr.theta.ex - pi/3 )
            
            invar.theta <- expr.theta.permu[ expr.theta.permu > -pi &
                expr.theta.permu <= pi ]
            invar.theta.ex <- c( invar.theta, 
                invar.theta + ifelse( invar.theta > 0, -2*pi, 2*pi ) )
            
            invar.theta.density <- density( invar.theta.ex, 
                adjust=theta.density.adjust, n=theta.density.n )
            invar.theta.density.sel <- invar.theta.density$x > -pi & 
                invar.theta.density$x <= pi
            
            # calculate Watson U2 statistic
            theta.watson.u2[ dim.group.idx ] <- 
                watson.u2( expr.theta, sample( invar.theta, 
                    length( expr.theta ), replace=TRUE ) )
        }
    }
    
    stdvec.result <- list( value = stdvec.offset, 
        stderr = stdvec.offset.stderr, 
        numerical.demand = stdvec.numerical.demand, 
        watson.u2 = theta.watson.u2, 
        h0.probe = h0.probe )
    
    stdvec.result
}

