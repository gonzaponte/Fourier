'''
    Module with useful tools.
    
    Author: G. Martinez Lema
    Date  : 17/10/2014
    
'''
from __future__ import division
import math
import numpy
import ROOT
import NEXT

class Array( numpy.ndarray ):
    '''
        Numpy array with some useful operations to apply on them.
    '''
    
    def __init__( self, arr ):
        '''
            Construct array from list of lists or copy array.
        '''
        print len(arr),len(arr[0])
        self.array = numpy.array( arr )
    
    def Sum( self, absvalue = abs ):
        '''
            Sum the content of the array. If absvalue = abs (default) the absolute value is taken. Use absvalue = None to avoid this behaviour.
        '''
        if absvalue is None:
            absvalue = lambda x: x
        return sum( sum( absvalue(self.array) ) )
    
    def Average( self ):
        '''
            Compute the weigthed average x,y-bin of the array.
        '''
        imean = jmean = 0.
        for i, row in enumerate( self.array ):
            for j, value in enumerate( row ):
                imean += i * value
                jmean += j * value
        
        total = self.Sum( None )
        return imean / total, jmean / total

    def Max( self ):
        '''
            Find greatest element.
        '''
        greatest = 0.
        i_greatest = j_greatest = 0
        for i, row in enumerate( self.array ):
            for j, value in enumerate( row ):
                if abs(value) > greatest:
                    greatest = abs(value)
                    i_greatest = i
                    j_greatest = j
        
        return i_greatest, j_greatest, greatest
    
    def CountNonZeros( self, eps = 1e-10 ):
        '''
            Count the number of elements different from 0 or eps.
        '''
        n = 0
        for row in self.array:
            for value in row:
                if value and abs(value) > eps:
                    n += 1
        return n
    
    def MakePSF( self ):
        '''
           Build the PSF graph.
        '''
        ib, jb = self.Average()
        qtot   = self.Sum()
        graph = ROOT.TGraph()
        n = 0
        for i,row in enumerate(self.array):
            for j,value in enumerate(row):
                if not value:
                    continue
                dr = math.sqrt( ( i - ib )**2 + ( j - jb )**2 ) * NEXT.TrackingPlane.Pitch
                q  = value / qtot
                graph.SetPoint( n, dr, q )
                n += 1
        return graph
    
    def FillHistogram( self, title = 'array', xlabel = 'x (mm)', ylabel = 'y (mm)', zlabel = 'Q (pes)' ):
        '''
            Fill an histogram with the data of the array.
        '''
        hxy = NEXT.TrackingPlane.GetHistogram( title, title )
        hxy.GetXaxis().SetTitle( xlabel ); hxy.GetXaxis().CenterTitle()
        hxy.GetYaxis().SetTitle( ylabel ); hxy.GetYaxis().CenterTitle()
        
        for i,row in enumerate(self.array):
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
        return Array( [ [ histogram.GetBinContent( i, j ) for j in range( 1, NEXT.TrackingPlane.Nbins + 1 ) ]
                                                          for i in range( 1, NEXT.TrackingPlane.Nbins + 1 ) ] )


if __name__ == '__main__':
    f = ROOT.TFile('Output/structures_0mm.root')
    h = f.Get('Uniform circle_xy;3')
    full = Array.FromHistogram(h)
    disc = NEXT.TrackingPlane.Discretize(full)
    psf = disc.MakePSF()
    psf.SetMarkerStyle(20)
    psf.SetMarkerSize(1)
#    fun = ROOT.TF1('fit','[0]/( 1 + [1]*x^2 + [2]*x^4 )^1.5')
    fun = ROOT.TF1('fit','[0]/( 1 + [1]*x^2 )^1.5')
#    fun.SetParameters( 1,1e-3, 1e-6 )
    fun.SetParameters( 1,1e-3 )
    fun.SetParLimits(0,0,10)
    fun.SetParLimits(1,0,10)
#    fun.SetParLimits(2,0,10)
#    fun.SetParLimits(3,0,10)
    psf.Fit( fun )
    psf.Draw('AP')
