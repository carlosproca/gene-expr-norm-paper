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


# Implements median condition-decomposition normalization


normalize.median.condec <- function( expression.data, expression.condition, 
    norm.probability=0.5, verbose=FALSE, p.value.graph=NULL )
{
    # identify samples and conditions to normalize
    normalize.sample <- 
        colnames( expression.data )[ ! is.na( expression.condition ) ]
    normalize.sample.condition <- 
        as.character( na.omit( expression.condition ) )
    normalize.condition <- unique( normalize.sample.condition )
    
    if ( length( normalize.condition ) == 0 )
        stop( "No condition to normalize" )
    else if ( min( table( normalize.sample.condition ) ) < 2 )
        stop( "There must be 2 or more samples in each condition" )
    
    # select probes for normalization
    # no missing values in any normalization sample
    expression.probe <- rownames( expression.data )
    normalize.probe.idx <- which( rowSums( 
        is.na( expression.data[ , normalize.sample ] ) ) == 0 )
    
    # normalize within conditions

    normalize.expr.data <- matrix( nrow = length( expression.probe ), 
        ncol = length( normalize.sample ) )
    rownames( normalize.expr.data ) <- expression.probe
    colnames( normalize.expr.data ) <- normalize.sample
    
    normalize.within.cond.offset <- rep( 0, length( normalize.sample ) )
    names( normalize.within.cond.offset ) <- normalize.sample
    
    normalize.within.cond.convergence <- vector( "list", 
        length( normalize.condition ) )
    names( normalize.within.cond.convergence ) <- normalize.condition
    
    condition.sample.idx <- split( 1 : length( expression.condition ), 
        expression.condition )
    
    for ( sample.idx in condition.sample.idx )
    {
        condition <- as.character( expression.condition[ sample.idx[ 1 ] ] )
        norm.sample.idx <- which( condition == normalize.sample.condition )
        
        within.cond.norm.result <- normalize.median.within.condition(
            expression.data[ , sample.idx ], normalize.probe.idx, condition, 
            norm.probability, verbose )
        
        normalize.expr.data[ , norm.sample.idx ] <- 
            within.cond.norm.result$data
        
        normalize.within.cond.offset[ norm.sample.idx ] <- 
            within.cond.norm.result$offset
    }
    
    # normalize between conditions
    
    normalize.between.cond.offset <- rep( 0, length( normalize.condition ) )
    names( normalize.between.cond.offset ) <- normalize.condition
    
    normalize.between.cond.h0.probe <- NULL
    normalize.between.cond.convergence <- NULL
    
    if ( length( normalize.condition ) > 1 )
    {
        between.cond.norm.result <- 
            normalize.median.between.condition( normalize.expr.data, 
                normalize.probe.idx, normalize.condition, 
                normalize.sample.condition, norm.probability, verbose, 
                p.value.graph )
        
        normalize.expr.data <- between.cond.norm.result$data
        
        normalize.between.cond.offset <- between.cond.norm.result$offset
        
        normalize.between.cond.h0.probe <- between.cond.norm.result$h0.probe
        
        normalize.between.cond.convergence <-
            between.cond.norm.result$convergence
    }
    
    normalize.offset <- normalize.within.cond.offset[ normalize.sample ] + 
        normalize.between.cond.offset[ normalize.sample.condition ]
    
    list( data = normalize.expr.data, offset = normalize.offset, 
        within.condition.offset = normalize.within.cond.offset, 
        between.condition.offset = normalize.between.cond.offset, 
        between.condition.h0.probe = normalize.between.cond.h0.probe, 
        between.condition.convergence = normalize.between.cond.convergence )
}


normalize.median.within.condition <- function( edata, norm.probe.idx, 
    condition, norm.prob, verbose )
{
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    # find median offset
    norm.median.offset <- apply( edata[ norm.probe.idx, ], 2, quantile, 
        probs=norm.prob )
    norm.median.offset <- norm.median.offset - mean( norm.median.offset )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.median.offset )
    
    list( data = edata, offset = norm.median.offset )
}


normalize.median.between.condition <- function( edata, norm.probe.idx, 
    norm.cond, norm.sample.cond, norm.prob, verbose, p.value.graph )
{
    # identify samples available per condition
    within.cond.n <- as.vector( table( norm.sample.cond )[ norm.cond ] )
    names( within.cond.n ) <- norm.cond
    
    # calculate balanced within-condition means
    bal.mean.n <- min( within.cond.n )
    
    expr.bal.mean.data <- sapply( norm.cond, function( cond ) 
        if ( within.cond.n[ cond ] == bal.mean.n )
            rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] )
        else  # within.cond.n[ cond ] > bal.mean.n
            rowMeans( t( apply( 
                t( edata[ norm.probe.idx, cond == norm.sample.cond ] ), 
                2, function( ed ) sample( ed, bal.mean.n ) ) ) )
    )
    
    # for F-statistic
    # calculate within-condition means
    expr.mean.data <- sapply( norm.cond, function( cond ) 
        rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] ) )
    
    # calculate within-condition variances
    within.cond.var <- sweep( sapply( norm.cond, function( cond ) 
        apply( edata[ norm.probe.idx, cond == norm.sample.cond ], 1, var ) ), 
        2, within.cond.n - 1, "*" )
    within.cond.var <- rowSums( within.cond.var ) / 
        ( sum( within.cond.n ) - length( norm.cond ) )
    
    # calculate normalization
    norm.median.result <- normalize.median.selection( expr.bal.mean.data, 
        expr.mean.data, within.cond.var, within.cond.n, norm.prob, verbose, 
        p.value.graph )
    
    # normalize each condition with its offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <-
            edata[ , cond == norm.sample.cond ] - 
            norm.median.result$offset[ cond ]
    
    list( data = edata, offset = norm.median.result$offset, 
        h0.probe = norm.median.result$h0.probe, 
        convergence = norm.median.result$convergence )
}


normalize.median.selection <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, verbose, p.value.graph )
{
    iter.max <- 100
    median.offset.accum.step.max <- 10
    median.offset.accum.threshold <- 0.1
    median.offset.single.threshold <- 0.01
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    last.median.offset <- vector( "numeric" )
    last.median.h0.probe <- vector( "character" )
    
    norm.median.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.median.offset.ratio <- vector( "numeric" )
    norm.median.offset.accum.step <- vector( "numeric" )
    norm.median.h0.probe.num <- vector( "numeric" )

    median.offset <- rep( 0, ncol( edata ) )
    
    if ( verbose )
        cat( paste0( "between.condition\n" ) )
    
    iter <- 0
    median.offset.accum.step <- 0
    median.offset.ratio <- 1
    
    while ( iter < iter.max && 
        median.offset.accum.step < median.offset.accum.step.max &&
        median.offset.ratio >= median.offset.single.threshold )
    {
        iter <- iter + 1
        
        # obtain next step of median offset     
        median.offset.step <- calculate.median.offset( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, norm.prob, 
            p.value.graph, iter )
        
        median.offset.delta <- median.offset.step$value
        median.h0.probe <- median.offset.step$h0.probe
        
        median.h0.probe.num <- length( median.h0.probe )
        
        # check errors
        if ( any( is.nan( median.offset.delta ) ) )
            stop( "NaN error in normalize.median.condec" )
        
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
    
    if ( median.offset.accum.step < median.offset.accum.step.max && 
            median.offset.ratio >= median.offset.single.threshold )
        stop( "no convergence in normalize.median.condec" )
    
    dimnames( norm.median.offset ) <- NULL
    
    norm.median.convergence <- list( offset = norm.median.offset, 
        offset.ratio = norm.median.offset.ratio, 
        offset.accum.step = norm.median.offset.accum.step, 
        h0.probe.num = norm.median.h0.probe.num )
    
    list( offset = last.median.offset, h0.probe = last.median.h0.probe, 
        convergence = norm.median.convergence )
}


calculate.median.offset <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, p.value.graph, iter )
{
    ks.test.alpha <- 1e-3
    
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
    
    # identify H0 probes with one-sided up Kolmogorov-Smirnov test
    ks.test.d <- sqrt( - log( ks.test.alpha ) / 2 )
    
    epv <- sort( expr.p.value )
    epv.n <- length( expr.p.value )
    
    epv.i <- 1
    ks.test.D.up <- 1:epv.n / epv.n - epv
    ks.test.reject <- any( ks.test.D.up > ks.test.d / sqrt( epv.n ) )
    
    while ( ks.test.reject )
    {
        epv.i <- epv.i + 1
        
        ks.test.D.up <- 
            ( epv.i : epv.n - epv.i + 1 ) / ( epv.n - epv.i + 1 ) - 
            ( epv[ epv.i : epv.n ] - epv[ epv.i - 1 ] ) / 
            ( 1 - epv[ epv.i - 1 ] )
        
        ks.test.reject <- 
            any( ks.test.D.up > ks.test.d / sqrt( epv.n - epv.i + 1 ) )
    }
    
    epv.h0.i <- epv.i
    epv.h0.n <- epv.n - epv.i + 1
    epv.h0.p <-  ifelse( epv.i == 1, 0, epv[ epv.i - 1 ] )
    epv.h0.q <- ( epv.i - 1 ) / epv.n
    
    h0.probe.idx <- order( expr.p.value, decreasing=TRUE )[ 1 : epv.h0.n ]
    h0.probe <- names( expr.p.value )[ h0.probe.idx ]
    
    pi0.est <- ( 1 - epv.h0.q ) / ( 1 - epv.h0.p )
    
    # plot graph of p-values
    if ( ! is.null( p.value.graph ) )
    {
        if ( p.value.graph != "" )
        {
            if ( ! file.exists( p.value.graph ) )
                dir.create( p.value.graph, recursive=TRUE )
            
            png.filename <- sprintf( "%s/pvalue_iter%02d.png", 
                p.value.graph, iter )
            
            png( png.filename, width=640, height=360 )
        }
        
        par.default <- par( no.readonly=TRUE )
        
        par( mfrow = c(1,2), pty="s", xpd=FALSE, 
            mar = c( 3.4, 2.2, 2.2, 2.2 ), oma = c( 0, 3.8, 0, 0 ) )
        
        epv.quant <- 1:epv.n / epv.n
        epv.h0.idx <- epv.h0.i : epv.n
        
        epv.x0 <- epv.h0.p
        epv.y0 <- epv.h0.q
        epv.yd <- ks.test.d * ( sqrt( epv.h0.n ) / epv.n )
        
        xylim <- list( c(0,0), c( epv.h0.p, epv.h0.q ) )
        
        for ( i in 1:2  )
        {
            plot( 0, type="n", 
                xlim = c( xylim[[i]][1], 1 ), ylim = c( xylim[[i]][2], 1 ), 
                xlab="p-value", ylab="", cex.axis=1.3, cex.lab=1.5 )
            
            segments( x0 = c( 0, epv[ - epv.h0.idx ] ), 
                y0 = c( 0, epv.quant[ - epv.h0.idx ] ), 
                x1 = c( epv[ - epv.h0.idx ], epv[ epv.h0.idx ][ 1 ] ), 
                lwd=2 )
            segments( x0 = epv[ epv.h0.idx ], y0 = epv.quant[ epv.h0.idx ], 
                x1 = c( epv[ epv.h0.idx ][ -1 ], 1 ), lwd=2, col="red" )
            points( epv.h0.p, epv.h0.q, pch=20, cex=2.5 )
            
            segments( 0, 1 - pi0.est, 1, 1, col="blue", lwd=1.5, lty=2 )
            segments( 0, 1 - pi0.est + epv.yd, 1, 1 + epv.yd, col="blue", 
                lwd=1.5, lty=3 )
            
            if ( i==1 )
                graph.title <- sprintf( "iter=%02d", iter )
            else
                graph.title <- substitute( paste( "#", H[0], "=", h0.n ), 
                    list( h0.n = sprintf( "%5d", epv.h0.n ) ) )
            
            title( main=graph.title, line = ifelse( i==1, 1.4, 1.75 ), 
                font.main=1, cex.main=1.7 )
        }
        
        mtext( "F( p-value )", side=2, line=1.3, outer=TRUE, cex=1.5 )
        
        if ( p.value.graph != "" )
            dev.off()
        else
            par( par.default )
    }
    
    # identify probes for normalization
    median.probe <- h0.probe
    
    # calculate offset
    median.offset <- apply( edata[ median.probe, ], 2, quantile, 
        probs=norm.prob )
    median.offset <- median.offset - mean( median.offset )
    
    list( value = median.offset, h0.probe = h0.probe )
}

