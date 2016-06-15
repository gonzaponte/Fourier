'''
    Tools to perform (2D-) Fourier transforms.
    
    @Author: G. Martinez Lema
    @Date  : 16/10/2014
'''

from __future__ import division
import numpy
import math
import ROOT
import NEXT
import Tools

### Aliases
FT  = numpy.fft. fft2
IFT = numpy.fft.ifft2

SHIFT = numpy.fft.fftshift

FTS  = lambda x: SHIFT(  FT( x ) )
IFTS = lambda x: SHIFT( IFT( x ) )

Discretize = NEXT.TrackingPlane.Discretize

def Gauss( x, y, sigmax, sigmay ):
    '''
        Compute the value at x,y of a 2D normalized gaussian with mean 0,0 and standard deviations sigmax and sigmay.
    '''
    return 1 / ( 2 * math.pi * sigmax * sigmay ) * math.exp( - 0.5 * (  ( x/sigmax )**2 +  ( y/sigmay )**2 ) )

def GaussMatrix( sigmax = 8., sigmay = 8. ):
    '''
       Compute a matrix filled with a normalized gaussian of mean = 0,0 and standard deviations sigmax and sigmay.
    '''
    return numpy.array( [ [ Gauss( x, y, sigmax, sigmay ) for y in NEXT.TrackingPlane.Bins ] for x in NEXT.TrackingPlane.Bins ] )

def PSF( x, y, x0 = 0., y0 = 0., z0 = 4.5 ):
    '''
        Compute the value of the normalized PSF with parameters x0,y0,z0 at x,y.
    '''
    return math.pow( ( x-x0 )**2 + ( y-y0 )**2 + z0**2 , -1.5 ) / ( 2 * math.pi )

def PSFMatrix( x0 = 0., y0 = 0., z0 = 4.5 ):
    '''
        Compute a matrix filled with the normalized PSF with parameters x0,y0,z0.
    '''
    return numpy.array( [ [ PSF( x, y, x0, y0, z0 ) for y in NEXT.TrackingPlane.Bins ] for x in NEXT.TrackingPlane.Bins ] )
        
def DeltaMatrix( x0 = 0., y0 = 0 ):
    '''
        Compute a matrix filled with a delta function at x0,y0.
    '''
    matrix = numpy.zeros( ( NEXT.TrackingPlane.Nbins, NEXT.TrackingPlane.Nbins ) )
    matrix[ int( x0 + .5 * NEXT.TrackingPlane.Nbins ) ][ int( x0 + .5 * NEXT.TrackingPlane.Nbins ) ] = 1.
    return matrix

def Shift( array ):
    half = int(NEXT.TrackingPlane.Range)
    for i in range(half):
        for j in range(half):
            si, sj = i + half, j + half
            array[i][j], array[si][sj] = array[si][sj], array[i][j]
    for i in range(half,2*half):
        for j in range(half):
            si, sj = i - half, j + half
            array[i][j], array[si][sj] = array[si][sj], array[i][j]
    
class Transformer:
    '''
        Contains algorithms for source computing.
    '''
    def __init__( self, signal, sigma = 8., true = None ):
        '''
            Construct with the signal and optionally, with the true (full) information.
        '''
        self.signal = signal
        self.sigma  = sigma
        self.true   = true
        self.Gauss  = FT( GaussMatrix( sigma, sigma ) ) if sigma else None
        self.PSFT   = FT( PSFMatrix() )
    
        self.Source()
        self.Image()
        self.signaldif = abs( Discretize(self.image) - self.signal      )
        self.truedif   = abs(            self.image  - self.trueimage   )
    
    def Source( self ):
        '''
            Compute source.
        '''
        if self.sigma: ### Convolute with a gaussian
            self.source = IFT( FT(self.signal) * self.Gauss / self.PSFT )
        else:     ### Do not convolute with a gaussian
            self.source = IFTS( FT(self.signal) / self.PSFT )

        if not self.true is None:
            self.truesource = IFT( FT(self.true) * FT( GaussMatrix( 1., 1. ) ) / self.PSFT )

    def Image( self ):
        '''
            Compute the image that would produce the reconstructed source.
        '''
        self.image = IFTS( FT( self.source ) *  self.PSFT )
        
        if not self.true is None:
            self.trueimage = IFTS( FT( self.truesource ) *  self.PSFT )

    def Chi2( self ):
        '''
            Compute the chi2 to compare the real signal with the expected one.
        '''
        measured = self.signal
        expected = self.image
        chi2     = 0.
        for i in NEXT.TrackingPlane.SiPMs_indices:
            for j in NEXT.TrackingPlane.SiPMs_indices:
                if expected[i][j]:
                    chi2 += ( expected[i][j] - measured[i][j] )**2 / expected[i][j]

        return chi2

    def Plot( self ):
        '''
            Plot all data.
        '''
        canvas = ROOT.TCanvas()

        if not self.true is None:
            canvas.Divide(2,3)
            self.histos = [ Tools.Arrays.FillHistogram(     self.signal     , 'signal'      ),
                            Tools.Arrays.FillHistogram(     self.true       , 'true'        ),
                            Tools.Arrays.FillHistogram( abs(self.source    ), 'source'      ),
                            Tools.Arrays.FillHistogram( abs(self.truesource), 'true source' ),
                            Tools.Arrays.FillHistogram(     self.signaldif  , 'signal diff' ),
                            Tools.Arrays.FillHistogram(     self.truedif    , 'true diff'   )]
        else:
            canvas.Divide(2,2)
            self.histos = [ Tools.Arrays.FillHistogram(     self.signal     , 'signal'      ),
                            Tools.Arrays.FillHistogram( abs(self.source    ), 'source'      ),
                            Tools.Arrays.FillHistogram( abs(self.image     ), 'from source' ),
                            Tools.Arrays.FillHistogram(     self.signaldif  , 'signal diff' )]
        
        for i,h in enumerate(self.histos):
            canvas.cd(i+1)
            if h:
                h.Draw('zcol')
    
        canvas.Update()
        return canvas



t = Transformer( Discretize(PSFMatrix()), 0., PSFMatrix() )
l = t.Plot()



