#
#  Utils.py
#
#  Description:
#  Useful functions.
#
#  Author:  Gonzalo
#  Date:   19/01/2014
#
#  Last update 05/02/2014 <- Gonzalo#




import os.path
import math
import ROOT
import time
import copy

# Sorts a list relative to a given number.
Order = lambda L,number: list( zip( *sorted( zip( map( lambda x: abs(x-number), L ), L ) ) )[1]) 



# Sorts a second list with the order of the first one
def OrderLike( L1, L2, order = 0 ):
    r = list( zip( *sorted( zip( L1, L2 ) ) )[1] )
    if order:
        return Reverse(r)
    return r

# Inverts the order of a list without destroying the original
Reverse = lambda L: L[::-1]



# Inverts zipping: (a,b) -> (b,a)
Izip  = lambda L: zip( Reverse( zip(*L) ) )



# Asymmetry factor between two values
Asymmetry = lambda x,y: ( x - y )/( x + y )



#Returns the index of the greatest/lowest value in a list
MaxIndex  = lambda x: x.index( max(x) )
MinIndex  = lambda x: x.index( min(x) )



# Return an appropiate name for the output file in order to not overwrite them.
# Be careful: use full paths.
def GetOutputName( name ):
    i=0
    while True:
        if not os.path.exists( name + str(i) + '.root' ):
            return name + str(i) + '.root'
        i += 1



# Returns cumulative of a list.
Cumulative = lambda L: [ sum( L[ :(i+1) ] ) for i in range( len(L) ) ]



# Returns a vector (matrix) of size rows (x cols) filled with zeros
def zeros(rows,cols=False):
    if cols:
        return [ zeros( cols ) for i in range(rows) ]
    return [ 0.0 for i in range(rows) ]



# Returns a vector (matrix) of size rows (x cols) filled with ones
def ones(rows,cols=False):
    if cols:
        return [ ones( cols ) for i in range(rows) ]
    return [ 1.0 for i in range(rows) ]



# Returns a N x N identity matrix
def identity( N ):
    M = zeros( N, N )
    for i in range( N ):
        M[i][i]=1.0
    return M



# Sign function
sign = lambda x: -1 if x<0 else 1



# Performs matrix multiplication
def MxM( M1, M2 ):
    f1 = len( M1 )
    c1 = len( M1[0] )
    c2 = len( M2[0] )
    
    M = zeros( f1, c2 )
    
    for i in range(f1):
        for j in range(c2):
            for k in range(c1):
                M[i][j]+= M1[i][k]*M2[k][j]
    return M



# Performs transposition of a matrix
def transpose( M ):
    return [ [ M[j][i] for j in range(len(M[0])) ] for i in range(len(M)) ]



# Matrix diagonalization
def Diagonalize( M, p=1e-4 ):
    def findmax(MM):
        maximum = 0.
        maxpos  = [0,0]
        for i in range(len(MM)):
            for j in range(len(MM[0])):
                if i==j:
                    continue
                if abs(MM[i][j]) > maximum:
                    maximum = abs(MM[i][j])
                    maxpos  = [i,j]
        return maxpos
    
    def check( MM ):
        for i in range( len(MM) ):
            for j in range( i+1,len(MM[0]) ):
                if abs( MM[i][j] ) > p :
                    return False
        return True

    D  = copy.deepcopy( M )
    V  = identity( len(M) )

    while True:
        if check(D):
            for i in range( len(D) ):
                D[i][i] = abs(D[i][i])
            break
        
        row,col = findmax(D)
        t = ( D[col][col] - D[row][row] )/( 2.*D[row][col] )
        t = sign( t )/( abs(t) + math.sqrt( t**2 + 1 ) )
        c = 1./math.sqrt( t**2 + 1 )
        s = c * t
        R = identity( len(M) )
        R[row][row] =  c
        R[col][col] =  c
        R[row][col] =  s
        R[col][row] = -s
        D = MxM( transpose( R ), MxM( D, R ) )
        V = MxM( V, R )
    
    return D,V



# Calculates mean with or without weights
def Mean( data, weights=None ):
    if not weights:
        return sum(data)/float(len(data))
    
    return sum(map(lambda x,y: x*y, data, weights))/sum(weights)



# Calculates variance with or without weights
def Variance( data, weights = None, mean = None ):
    return Covariance( data, data, weights, mean, mean )



# Calculates RMS with or without weights
def RMS( data, weights = None, mean = None ):
    return math.sqrt( Variance( data, weights, mean ) )



# Calculates weights for a list of uncertainties
def Weights( uncertainties ):
    return map( lambda x: x**(-2.), uncertainties )



# Calculates the median for a given data.
def Median( data ):
    if isinstance( data, list ):
        if len(data)%2:
            return sorted(data)[len(data)/2]
        return 0.5*( sorted(data)[len(data)/2-1] + sorted(data)[len(data)/2] )
    
    data = sorted(data.items())
    acum = Acumulate( data )
    tot  = sum(zip(*data)[1])
    if tot/2 in acum:
        return d[acum.index(tot/2)]
    
    f1 = filter( lambda x: x>tot/2, acum )[0]
    indice = acum.index(f1)
    f0 = acum[indice-1]
    a1 = d[indice]
    a0 = d[indice-1]
    return a0 + (.5*tot - f0)/(f1-f0)*(a1-a0)



# Returns the mode of a distribution of data.
def Mode( data ):    
    if isinstance( data, ( list, tuple ) ):
        return Mode( Frecs( data ) )
    
    value = -1
    maximum = -float('inf')
    for i,j in data.items():
        if j>maximum:
            maximum = j
            value = i
    return value


# Calculates the covariance between a set of data
def Covariance( x, y, weights = None, xmean = None, ymean = None ):
    if len(x) == len(y) == 1:
        return 0.
    
    if not xmean or not ymean:
        xmean = Mean( x, weights )
        ymean = Mean( y, weights )

    if not weights:
        return sum( map( lambda x,y: (x-xmean)*(y-ymean), x, y ) )/float( len(x) )

    normfactor = sum( weights )**2 - sum ( map( lambda x: x**2, weights ) )
    normfactor = sum( weights )/normfactor

    return normfactor * sum( map( lambda x,y,w: w*(x-xmean)*(y-ymean), x, y, weights ) )



# Returns the correlation parameter between 2 sets of data.
def Correlation( x, y, weights = None, xmean = None, ymean = None ):
    if len(x) == len(y) < 2:
        return 1. if x==y else 0.
    
    return Covariance( x, y, weights, xmean, ymean ) / ( RMS( x, weights, xmean ) * RMS( y, weights, ymean ) )


def ElapsedTime( t0 = 0 ):
    dt  = time.time() - t0
    h   = dt // 3600
    dt -= h * 3600
    m   = dt // 60
    dt -= m * 60
    s   = dt
    
    return '{0} h + {1} min + {2} s'.format( h, m, s )

def FillHisto( matrix, title ):
    '''Fills a histogram with a matrix.'''
    histo = ROOT.TH2F( title, title, 151, -75.5, 75.5, 151, -75.5, 75.5 )
    histo.GetXaxis().SetTitle('x (mm)')
    histo.GetYaxis().SetTitle('y (mm)')
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    for i in range( 151 ):
        for j in range( 151 ):
            val = abs(matrix[i][j])
            histo.SetBinContent( histo.GetBin( i+1, j+1 ), val )
            histo.SetBinError  ( histo.GetBin( i+1, j+1 ), math.sqrt(val) )
    return histo

def Resolution( sigma, mean ):
    return 2.35 * sigma / float(mean) * 100, 2.35 * sigma / float(mean * 48102)**0.5 * 100

def Scale2D( M, F ):
    '''Scales a 2D matrix by a factor F.'''
    return [ [ F * M[i][j] for j in range(len(M[0])) ] for i in range(len(M)) ]

def Sum2D( M, absval = False ):
    '''Sums up a 2D matrix by a factor F.'''
    return sum( map( sum, Abs2D( M ) ) ) if absval else sum( map( sum, M ) )

def Abs2D( M ):
    '''Absolute value of a 2D matrix.'''
    return [ [ abs( M[i][j] ) for j in range(len(M[0])) ] for i in range(len(M)) ]

def Multiply2D( M1, M2 ):
    '''Scales a 2D matrix by a factor F.'''
    return [ [ M1[i][j] * M2[i][j] for j in range(len(M1[0])) ] for i in range(len(M2)) ]

def Max2D( M ):
    '''Finds the greatest element of a 2D matrix.'''
    
    greatest = 0.
    pos      = 0, 0
    for x,col in enumerate(M):
        for y,value in enumerate(col):
            value = abs(value)
            if value > greatest:
                greatest = value
                pos = x,y
    
    return pos[0]-75,pos[1]-75,greatest

def Average2D( M ):
    '''Averages the data of a 2D matrix.'''
    
    xm, ym, I = 0., 0., 0.
    
    for x,col in enumerate(M):
        x -= 75
        for y,value in enumerate(col):
            y -= 75
            value = abs(value)
            
            xm += x * value
            ym += y * value
            I  +=     value
    
    return xm/I, ym/I

def ComputeFWHM( H ):
    binmax  = H.GetMaximumBin()
    fwhmval = H.GetBinContent( binmax ) * .5
    binlow  = binmax - 1
    binupp  = binmax + 1
    
    while H.GetBinContent( binlow ) > fwhmval:
        binlow -= 1
    while H.GetBinContent( binupp ) > fwhmval:
        binupp -= 1

    return H.GetBinCenter( binupp ) - H.GetBinCenter( binlow )

def Angle( x, y ):
    phi  = math.atan( y/x )
    if x < 0:
        phi += math.pi if x < 0 else 2 * math.pi if y < 0 else 0.
    return phi

def Corona( r, r0 ):
    return 1 + r // r0

def PhiSector( phi, phi0 ):
    return 1 + phi // phi0

