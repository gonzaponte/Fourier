'''
    Tools to perform (2D-) Fourier transforms.
    
    Author: G. Martinez Lema
    Date  : 16/10/2014
'''

from __future__ import division
from FourierToy import TrackingPlane
import numpy
import math
import ROOT


### Aliases
FT  = numpy.fft. fft2
IFT = numpy.fft.ifft2

SHIFT = numpy.fft.fftshift

FTS  = lambda x: SHIFT(  FT( x ) )
IFTS = lambda x: SHIFT( IFT( x ) )


class ArrayOperations:
    '''
        Some useful operations to apply on 2D-arrays.
    '''
    
    @staticmethod
    def Sum( array, absvalue = abs ):
        '''
            Sum the content of the array. If absvalue = abs (default) the absolute value is taken. Use absvalue = None to avoid this behaviour.
        '''
        if absvalue is None:
            absvalue = lambda x: x
        return sum( sum( absvalue(array) ) )

    @staticmethod
    def Average( array ):
        '''
            Compute the weigthed average x,y-bin of the array.
        '''
        imean = jmean = 0.
        for i, row in enumerate( array ):
            for j, value in enumerate( row ):
                imean += i * value
                jmean += j * value

        total = ArrayOperations.Sum( array, None )
        return imean / total, jmean / total

    @staticmethod
    def Max( array ):
        '''
            Find greatest element.
        '''
        greatest = 0.
        i_greatest = j_greatest = 0
        for i, row in enumerate( array ):
            for j, value in enumerate( row ):
                if abs(value) > greatest:
                    greatest = abs(value)
                    i_greatest = i
                    j_greatest = j
    
        return i_greatest, j_greatest, greatest

    @staticmethod
    def CountNonZeros( array, eps = 1e-10 ):
        n = 0
        for row in array:
            for value in row:
                if value and abs(value) > eps:
                    n += 1
        return n

    @staticmethod
    def GetHistogram( array, title = 'array', xlabel = 'x (mm)', ylabel = 'y (mm)', zlabel = 'Q (pes)' ):
        hxy = ROOT.TH2D( title, title, TrackingPlane.Nbins, -TrackingPlane.Range, TrackingPlane.Range,
                                       TrackingPlane.Nbins, -TrackingPlane.Range, TrackingPlane.Range )
        
        hxy.GetXaxis().SetTitle( xlabel ); hxy.GetXaxis().CenterTitle()
        hxy.GetYaxis().SetTitle( ylabel ); hxy.GetYaxis().CenterTitle()
        
        for i,row in enumerate(array):
            for j,value in enumerate(row):
                bin = hxy.GetBin( i+1, j+1 )
                hxy.SetBinContent( bin , value )
                hxy.SetBinError  ( bin , math.sqrt( value ) )

        return hxy

    @staticmethod
    def FromHistogram( histogram ):
        '''
            Construct array with the histogram data.
        '''
        return np.array( [ [ histogram.GetBinContent( i, j ) for j in range( 1, TrackingPlane.Nbins + 1 ) ]
                                                             for i in range( 1, TrackingPlane.Nbins + 1 ) ] )

def Gauss( x, y, sigmax, sigmay ):
    '''
        Compute the value at x,y of a 2D normalized gaussian with mean 0,0 and standard deviations sigmax and sigmay.
    '''
    return 1 / ( 2 * math.pi * sigmax * sigmay ) * math.exp( - 0.5 * (  ( x/sigmax )**2 +  ( y/sigmay )**2 ) )

def GaussMatrix( sigmax = 8., sigmay = 8. ):
    '''
       Compute a matrix filled with a normalized gaussian of mean = 0,0 and standard deviations sigmax and sigmay.
    '''
    return numpy.array( [ [ Gauss( x, y, sigmax, sigmay ) for y in TrackingPlane.Pitch ] for x in TrackingPlane.Pitch ] )

def PSF( x, y, x0 = 0., y0 = 0., z0 = 4.5 ):
    '''
        Compute the value of the normalized PSF with parameters x0,y0,z0 at x,y.
    '''
    return math.pow( ( x-x0 )**2 + ( y-y0 )**2 + z0**2 , -1.5 ) / ( 2 * math.pi )

def PSFMatrix( x0 = 0., y0 = 0., z0 = 4.5 ):
    '''
        Compute a matrix filled with the normalized PSF with parameters x0,y0,z0.
    '''
    return numpy.array( [ [ PSF( x, y, x0, y0, z0 ) for y in TrackingPlane.Pitch ] for x in TrackingPlane.Pitch ] )
        
def DeltaMatrix( x0 = 0., y0 = 0 ):
    '''
        Compute a matrix filled with a delta function at x0,y0.
    '''
    matrix = numpy.zeros( ( TrackingPlane.Nbins, TrackingPlane.Nbins ) )
    matrix[ int( x0 + .5 * TrackingPlane.Nbins ) ][ int( x0 + .5 * TrackingPlane.Nbins ) ] = 1.
    return matrix

def ComputeSource( I, sigma = 8. ):
    '''
        Compute the source that produced the signal I in the tracking plane.
    '''
    return IFT( FT(I) * FT( GaussMatrix( sigma, sigma ) ) / FT( PSFMatrix() ) ) )

def ComputeImage( S ):
    '''
        Compute the source that produced the signal I in the tracking plane.
    '''
    return SHIFT( IFT( FT( S ) *  FT( PSFMatrix() ) ) )

def ComputeChi2( signal, source ):
    '''
        Compute the chi2 to compare the real signal with the expected one.
    '''
    expected = ComputeImage( source )
    chi2 = 0.
    for i in range(len(expected)):
        for j in range(len(expected)):
            if expected[i][j]:
                chi2 += ( expected[i][j] - signal[i][j] )**2 / expected[i][j]

    return chi2






