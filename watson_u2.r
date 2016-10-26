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


# Implements Watson U2 statistic

watson.u2 <- function( x, y )
{
    # see section 6.5 of Durbin, Distribution Theory for Tests Based on the 
    # Sample Distribution Function, SIAM, Philadelphia (1973)
    
    n <- length( x )
    m <- length( y )
    
    r <- c( sort( x ), sort( y ) )
    r.rank <- rank( r, ties.method="average" )
    
    z <- ( r.rank[ 1:n ] - 1:n ) / m - ( 1:n - 1/2 ) / n 
    
    ( m / (n+m) ) * sum( ( z - mean( z ) )^2 ) + 
        ( m*(m+2*n) ) / ( 12*n*m*(n+m) ) 
}

