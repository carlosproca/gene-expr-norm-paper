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


# Adapts standard-vector condition-decompositon normalization for study across 
# probes


normalize.standard.vector.probe <- function( edata, edata.fstat, 
    within.cond.var, within.cond.n, condition, h0.probe.random, h0.probe.num, 
    verbose, vector.graph )
{
    iter.max <- 100
    stdvec.offset.accum.step.max <- 10
    stdvec.offset.accum.threshold <- 0.5    # 0.1 -> 0.5 increase needed for 
    stdvec.offset.single.threshold <- 0.01  # convergence with very small 
    vector.graph.probe.num <- 5000          # h0.probe.num
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    last.stdvec.offset <- vector( "numeric" )
    last.stdvec.h0.probe <- vector( "character" )
    
    norm.stdvec.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.stdvec.offset.sd <- vector( "numeric" )
    norm.stdvec.offset.stderr <- vector( "numeric" )
    norm.stdvec.offset.delta.sd <- vector( "numeric" )
    norm.stdvec.offset.accum.step <- vector( "numeric" )
    norm.stdvec.numerical.demand <- vector( "numeric" )
    norm.stdvec.watson.u2 <- matrix( nrow=0, 
        ncol = ( ncol( edata ) - 1 ) %/% 3 + 1 )
    norm.stdvec.h0.probe.num <- vector( "numeric" )
    
    stdvec.offset <- rep( 0, ncol( edata ) )
    
    if ( h0.probe.random )
        h0.probe <- sample( rownames( edata.fstat ), h0.probe.num )
    else
        h0.probe <- NULL
    
    vector.graph.probe <- NULL
    
    if ( ! is.null( vector.graph ) )
    {
        # select a sample of probes for plotting standardized sample vectors
        edata.probe <- rownames( edata )
        edata.probe.num <- length( edata.probe )
        
        if ( edata.probe.num <= vector.graph.probe.num )
            vector.graph.probe <- edata.probe
        else
            vector.graph.probe <- sample( edata.probe, vector.graph.probe.num )
    }
    
    iter <- 0
    stdvec.offset.accum.step <- 0
    stdvec.offset.ratio <- 1
    
    while ( iter < iter.max && 
        stdvec.offset.accum.step < stdvec.offset.accum.step.max &&
        stdvec.offset.ratio >= stdvec.offset.single.threshold )
    {
        iter <- iter + 1
        
        # obtain next step of standard vector offset     
        stdvec.offset.step <- calculate.stdvec.offset.probe( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, h0.probe, 
            h0.probe.num, vector.graph, vector.graph.probe, condition, iter )
        
        stdvec.offset.delta <- stdvec.offset.step$value
        stdvec.offset.stderr <- stdvec.offset.step$stderr
        stdvec.numerical.demand <- stdvec.offset.step$numerical.demand
        stdvec.watson.u2 <- stdvec.offset.step$watson.u2
        stdvec.h0.probe <- stdvec.offset.step$h0.probe
        
        stdvec.h0.probe.num <- NULL
        if ( ! is.null( stdvec.h0.probe ) )
            stdvec.h0.probe.num <- length( stdvec.h0.probe )
        
        # check errors
        if ( any( is.nan( stdvec.offset.delta ) ) || 
                is.nan( stdvec.offset.stderr ) )
            stop( "NaN error in normalize.stdvec.condec.probe" )
        
        if ( stdvec.numerical.demand < .Machine$double.eps * 10^3 )
            stop( "numerical error in normalize.stdvec.condec.probe" )
        
        # update total standard vector offset
        stdvec.offset <- stdvec.offset + stdvec.offset.delta
        stdvec.offset <- stdvec.offset - mean( stdvec.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, stdvec.offset )
        if ( ! is.null( edata.fstat ) )
            last.norm.data.fstat <- sweep( edata.fstat, 2, stdvec.offset )
        
        # check convergence
        stdvec.offset.sd <- sd( stdvec.offset )
        stdvec.offset.delta.sd <- sd( stdvec.offset.delta )
        
        stdvec.offset.ratio <- stdvec.offset.delta.sd / stdvec.offset.stderr
        
        if ( stdvec.offset.ratio < stdvec.offset.accum.threshold )
            stdvec.offset.accum.step <- stdvec.offset.accum.step + 1
        else
            stdvec.offset.accum.step <- 0
        
        # store last results
        last.stdvec.offset <- stdvec.offset
        last.stdvec.h0.probe <- stdvec.h0.probe
        
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
        if ( ! is.null( stdvec.watson.u2 ) )
            norm.stdvec.watson.u2 <- rbind( norm.stdvec.watson.u2, 
                stdvec.watson.u2 )
        if ( ! is.null( stdvec.h0.probe.num ) )
            norm.stdvec.h0.probe.num <- c( norm.stdvec.h0.probe.num, 
                stdvec.h0.probe.num )
        
        if ( verbose )
        {
            stdvec.offset.stderr.ratio <- 
                stdvec.offset.stderr / stdvec.offset.sd
            stdvec.offset.delta.sd.ratio <- 
                stdvec.offset.delta.sd / stdvec.offset.sd
            
            if ( ! is.null( stdvec.watson.u2 ) )
                stdvec.watson.u2.char <- paste0( signif( stdvec.watson.u2, 6 ), 
                    collapse=" " )
            else
                stdvec.watson.u2.char <- ""
            
            cat( sprintf( "  %2d %g %g %g %02d %g [%s]", iter, 
                stdvec.offset.sd, stdvec.offset.stderr.ratio, 
                stdvec.offset.delta.sd.ratio, stdvec.offset.accum.step, 
                stdvec.numerical.demand, stdvec.watson.u2.char ), 
                stdvec.h0.probe.num, "\n" )
        }
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( iter == iter.max )
        stop( "no convergence in normalize.stdvec.condec.probe" )
    
    dimnames( norm.stdvec.offset ) <- NULL
    dimnames( norm.stdvec.watson.u2 ) <- NULL

    norm.stdvec.convergence <- list( offset = norm.stdvec.offset, 
        offset.sd = norm.stdvec.offset.sd, 
        offset.stderr = norm.stdvec.offset.stderr, 
        offset.delta.sd = norm.stdvec.offset.delta.sd, 
        offset.accum.step = norm.stdvec.offset.accum.step, 
        numerical.demand = norm.stdvec.numerical.demand, 
        watson.u2 = norm.stdvec.watson.u2 )
    
    if ( length( norm.stdvec.h0.probe.num ) > 0 )
        norm.stdvec.convergence$h0.probe.num <- norm.stdvec.h0.probe.num
    
    if ( is.null( last.stdvec.h0.probe ) ) {
        list( offset = last.stdvec.offset, 
            convergence = norm.stdvec.convergence )
    } else {
        list( offset = last.stdvec.offset, h0.probe = last.stdvec.h0.probe, 
            convergence = norm.stdvec.convergence )
    }
}


calculate.stdvec.offset.probe <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, h0.probe, h0.probe.num, vector.graph, vector.graph.probe, 
    condition, iter )
{
    stdvec.trim = 0.01
    
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
    expr.var <- apply( edata.fstat[ h0.probe, ], 1, var )
    expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
    
    stdvec.probe <- h0.probe[ 
        expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
        expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ]
    
    stdvec.probe <- intersect( stdvec.probe, rownames( edata ) )
    
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
    
    # calculate and plot density distribution of standard vector angles
    
    dimension.num <- ncol( expr.scaled )
    
    if ( dimension.num < 3 )
    {
        theta.watson.u2 = NULL
    }
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
            
            if ( ! is.null( vector.graph ) && require( plotrix, quietly=TRUE ) )
            {
                if ( vector.graph != "" )
                {
                    if ( ! file.exists( vector.graph ) )
                        dir.create( vector.graph, recursive=TRUE )
                    
                    png.filename <- sprintf( 
                        "%s/stdvec_%s%s_iter%02d.png", vector.graph, condition, 
                        ifelse( dimension.group.num > 1, 
                            sprintf( "_dg%02d", dim.group.idx ), "" ), 
                        iter )
                    
                    png( png.filename, width=640, height=360 )
                }
                
                par.default <- par( no.readonly=TRUE )
                
                par( mfrow = c(1,2), pty="s", xpd=FALSE, 
                    mar = c( 1.2, 1.2, 2.2, 1.2 ), oma = c( 0, 0, 0, 0 ) )
                
                # select offset values for condition group
                stdvec.offset.uv <- stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ] %*% uv
                stdvec.offset.u <- stdvec.offset.uv[ , 1 ]
                stdvec.offset.v <- stdvec.offset.uv[ , 2 ]

                # plot a sample of standardized sample vectors
                uv.lim <- c( -1.5, 1.5 )
                plot( 0, type="n", xlim=uv.lim, ylim=uv.lim, axes=FALSE, 
                    ann=FALSE, frame.plot=FALSE, asp=1 )
                
                expr.uv.probe <- names( expr.u )
                if ( length( expr.uv.probe ) > length( vector.graph.probe ) ) {
                    expr.uv.probe <- intersect( expr.uv.probe, 
                        vector.graph.probe )
                }
                expr.uv.probe.num <- length( expr.uv.probe )
                
                expr.uv.factor <- 1.05
                expr.uv.color <- gray( 0.3 )
                expr.uv.width <- ifelse( expr.uv.probe.num > 1000, 0.1, 0.2 )
                segments( 0, 0, expr.uv.factor * expr.u[ expr.uv.probe ], 
                    expr.uv.factor * expr.v[ expr.uv.probe ], 
                    col=expr.uv.color, lwd=expr.uv.width )
                
                grid.pos.x <- c( 0, -sqrt(3/4), sqrt(3/4) )
                grid.pos.y <- c( 1, -1/2, -1/2 )
                
                grid.length <- expr.uv.factor * sqrt( 2 )
                segments( 0, 0, grid.length * grid.pos.x, 
                    grid.length * grid.pos.y, lwd=2 )
                
                stdvec.offset.uv.factor <- 10
                stdvec.offset.uv.color <- "red"
                segments( 0, 0, stdvec.offset.uv.factor * stdvec.offset.u, 
                    stdvec.offset.uv.factor * stdvec.offset.v, lwd=2, 
                    col=stdvec.offset.uv.color )
                
                par( xpd=TRUE )
                
                grid.label.length <- 1.63
                grid.labels <- c( "s1", "s2", "s3" )
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels )
                
                par( xpd=FALSE )
                
                stdvec.offset.mag <- sqrt( sum( stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ]^2 ) )
                mtext( paste0( "||offset|| = ", 
                    signif( stdvec.offset.mag, 4 ) ), side=1 )
                
                # plot polar distributions of standard vector angles
                polar.expr.theta <- expr.theta.density$x[ 
                    expr.theta.density.sel ]
                polar.expr.rho <- 2 * expr.theta.density$y[ 
                    expr.theta.density.sel ]

                polar.invar.theta <- invar.theta.density$x[ 
                    invar.theta.density.sel ]
                polar.invar.rho <- 2 * invar.theta.density$y[ 
                    invar.theta.density.sel ]
                
                rho.grid <- seq( 0, 3/(2*pi), length.out=4 )
                rho.labels <- c( "", "", expression(1/pi), "" )
                theta.labels <- c( "", "", "" )
                
                radial.plot( polar.expr.rho, polar.expr.theta - pi/2, 
                    start=pi/2, rp.type="p", radial.lim = rho.grid, 
                    radial.labels = rho.labels, 
                    show.grid.labels = length( theta.labels ), 
                    labels = theta.labels, mar = par( "mar" ) )
                
                radial.plot( polar.invar.rho, polar.invar.theta - pi/2, 
                    start=pi/2, rp.type="p", radial.lim = rho.grid, lty=2, 
                    line.col="blue", add=TRUE )
                
                grid.label.length <- 0.522
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels )
                
                mtext( substitute( paste( "Watson U"^"2", " = ", wu2 ), 
                    list( wu2 = sprintf( "%.3e", 
                        theta.watson.u2[ dim.group.idx ] ) ) ), 
                    side=1, line=0.2, cex=1.5 )
                
                graph.title <- sprintf( "%s%s    -    iter=%02d", condition, 
                    ifelse( dimension.group.num > 1, 
                        sprintf( ":%02d", dim.group.idx ), "" ), 
                    iter )
                title( main=graph.title, outer=TRUE, line=-1.8, font.main=1, 
                    cex.main=1.7 )
                
                if ( vector.graph != "" )
                    dev.off()
                else
                    par( par.default )
            }
        }
    }
    
    list( value = stdvec.offset, stderr = stdvec.offset.stderr, 
        numerical.demand = stdvec.numerical.demand, 
        watson.u2 = theta.watson.u2, 
        h0.probe = h0.probe )
}

