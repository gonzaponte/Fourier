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

ROOT.gStyle.SetOptStat('')

def Shift( array ):
    half = int(NEXT.TrackingPlane.Range)
    for i in range(half):
        si = i + half
        for j in range(half):
            sj = j + half
            array[i][j], array[si][sj] = array[si][sj], array[i][j]
    for i in range(half,2*half):
        si = i - half
        for j in range(half):
            sj = j + half
            array[i][j], array[si][sj] = array[si][sj], array[i][j]
    return array


### Aliases
FT  = numpy.fft. fft2
IFT = numpy.fft.ifft2

SHIFT = numpy.fft.fftshift
SHIFT = Shift

FTS  = lambda x: SHIFT(  FT( x ) )
IFTS = lambda x: SHIFT( IFT( x ) )

Integral   = Tools.Arrays.Sum
Discretize = NEXT.TrackingPlane.Discretize
FromHistogram = Tools.Arrays.FromHistogram
FillHistogram = Tools.Arrays.FillHistogram

def plot( x ):
    c = ROOT.TCanvas()
    h = FillHistogram(abs(x))
    h.Draw('zcol')
    c.Update()
    return c,h


### Distributions

def PSF_data ( p0, p1, p2, p3):
    return lambda x,y,x0,y0,z0,N:  N * ( p0/1e3 * ( 1 + p1/1e3 * ( (x-x0)**2 + (y-y0)**2 ) + p2/1e6 * ( (x-x0)**2 + (y-y0)**2 )**2 )**-1.5 + p3/1e3 ) / 1.2009487822

PSF_Mip  = PSF_data( 71.69, 2.380, 2.406, 1.0000 )

def Gauss( x, y, sigmax, sigmay ):
    '''
        Compute the value at x,y of a 2D normalized gaussian with mean 0,0 and standard deviations sigmax and sigmay.
    '''
    return 1 / ( 2 * math.pi * sigmax * sigmay ) * math.exp( - 0.5 * (  ( x/sigmax )**2 +  ( y/sigmay )**2 ) )

def GaussMatrix( sigmax = 8., sigmay = 8. ):
    '''
       Compute a matrix filled with a normalized gaussian of mean = 0,0 and standard deviations sigmax and sigmay.
    '''
    return numpy.array( [ [ Gauss( x, y, sigmax, sigmay ) for y in NEXT.TrackingPlane.BinsS ] for x in NEXT.TrackingPlane.BinsS ] )

def PSF( x, y, x0 = 0., y0 = 0., z0 = 4.5, N = 1e8 ):
    '''
        Compute the value of the normalized PSF with parameters x0,y0,z0 at x,y.
    '''
    return math.pow( ( x-x0 )**2 + ( y-y0 )**2 + z0**2 , -1.5 ) * N * z0 / ( 2 * math.pi )

def PSFMatrix( x0 = 0., y0 = 0., z0 = 4.5, N = 1e9, psf = PSF ):
    '''
        Compute a matrix filled with the normalized PSF with parameters x0,y0,z0.
    '''
    return numpy.array( [ [ psf( x, y, x0, y0, z0, N ) for y in NEXT.TrackingPlane.BinsS ] for x in NEXT.TrackingPlane.BinsS ] )
        
def DeltaMatrix( x0 = 0., y0 = 0 ):
    '''
        Compute a matrix filled with a delta function at x0,y0.
    '''
    matrix = numpy.zeros( ( NEXT.TrackingPlane.NbinsS, NEXT.TrackingPlane.NbinsS ) )
    matrix[ int( x0 + .5 * NEXT.TrackingPlane.NbinsS ) ][ int( x0 + .5 * NEXT.TrackingPlane.NbinsS ) ] = 1.
    return matrix

def Source( signal, sigma = 8.):
    '''
        Compute source.
    '''
    if sigma: ### Convolute with a gaussian
        source = IFT( FT(signal) * FT( GaussMatrix( sigma, sigma ) ) / psfT )
    else:     ### Do not convolute with a gaussian
        source = IFTS( FT(signal) / psfT )
    return source

def Image( source ):
    '''
        Compute the image that would produce the reconstructed source.
    '''
    image = IFTS( FT( source ) * psfT )
    dimage = Discretize(image)
    return dimage * Integral( image ) / Integral( dimage )

def Residuals( image, guess, nsigmas = 2. ):
    new = numpy.array( image )
    guess = abs( guess )
    for i in NEXT.TrackingPlane.SiPMsS_indices:
        for j in NEXT.TrackingPlane.SiPMsS_indices:
            new[i][j] -= guess[i][j]
            if new[i][j] < nsigmas * math.sqrt(image[i][j]):
                new[i][j]  = 0.
    return new

def Filter( source, cut = 0.99 ):
    '''
        Cut reconstructed source to cut level.
        '''
    newsource = numpy.array(abs(source))
    integral = Tools.Arrays.Sum( newsource )
    maxval = Tools.Arrays.Max( newsource )[2]
    for i in NEXT.TrackingPlane.BinsS:
        for j in NEXT.TrackingPlane.BinsS:
            if newsource[i][j]/maxval < cut:
                newsource[i][j] = 0.
    return newsource * integral / Tools.Arrays.Sum( newsource )

psf  = PSFMatrix()
#psf  = PSFMatrix(psf = PSF_Mip)
psfT = FT( psf )

class Transformer:
    '''
        Contains algorithms for source computing.
    '''
    def __init__( self, signal, true = None ):
        '''
            Construct with the signal and optionally, with the true (full) information.
        '''
        self.signal = signal
        self.true   = true
        self.Gauss  = FT( GaussMatrix( 8., 8. ) )
    
    def Compute( self, sigma = 8., cut = 0.99, nsigmas = 2. ):
        '''
            Compute everything.
        '''
        if not self.signal.any():
            self.sources = [ self.signal ]
            self.fsource = [ self.signal ]
            return self.signal
        self.sources   = [ Source( self.signal, sigma ) ]
        self.filtered  = [ Filter( self.sources[-1], cut ) ]
        self.fsource   = self.filtered[-1]
        self.images    = [ self.signal, Image( self.fsource ) ]
        self.residuals = [ Residuals( self.signal, self.images[-1], nsigmas ) ]
        
        while Tools.Arrays.Sum( self.residuals[-1] ):
#            a = plot( self.residuals[-1] )
#            raw_input()
            self.sources.append( Source( self.residuals[-1] ) )
#            a = plot( self.sources[-1] )
#            raw_input()
            self.filtered.append( Filter( self.sources[-1] ) )
            self.fsource = self.fsource + self.filtered[-1]
            self.images.append( Image( self.fsource ) )
            self.residuals.append( Residuals( self.signal, self.images[-1] ) )
#            a = plot(self.residuals[-1])
#            raw_input()
        return self.fsource

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

def Shrink( large ):
    short = numpy.ndarray( (150,150) )
    for i in range(NEXT.TrackingPlane.NbinsS):
        for j in range(NEXT.TrackingPlane.NbinsS):
            short[i][j] = large[i][j]
    return short

if __name__ == '__main__':
#simus  = 'Output/structures_0mm.root'
#simus  = ROOT.TFile(simus)
#
#pointtrue  = FromHistogram( simus.Get('Pointlike_xy;1') )
#circletrue = FromHistogram( simus.Get('Uniform circle_xy;1') )
#linetrue   = FromHistogram( simus.Get('Straight line_xy;1') )
#
#pointtrue = Shrink( pointtrue )
#circletrue = Shrink( circletrue )
#linetrue   = Shrink( linetrue )
#
#point  = Discretize( pointtrue )
#circle = Discretize( circletrue )
#line   = Discretize( linetrue )
#
#
#t = Transformer( circle )
#t.Compute()
#a = plot( t.fsource )

#i0  = Discretize( PSFMatrix() )
#i0  = circle
#s0  = Source( i0 )
#s1  = Filter( s0, 0.95 )
#i1  = Image( s1 )
#i2  = Discretize( i1 )
#i2 *= Integral( i1 ) / Integral( i2 )
#r   = Residuals( i0, i2, 2 )

    alpha = .0001 * 1
    s = (1-alpha) * PSFMatrix(-30) + alpha * PSFMatrix(30)
    s = Discretize(s)
    t = Transformer(s)
    f = t.Compute(8.,0.99)
    a = plot(f)

