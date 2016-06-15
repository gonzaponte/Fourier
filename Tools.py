'''
    Module with useful tools.
    
    Author: G. Martinez Lema
    Date  : 17/10/2014
    
'''

import numpy
import ROOT
import NEXT

class Arrays:
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
        
        total = Arrays.Sum( array, None )
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
        '''
            Count the number of elements different from 0.
        '''
        n = 0
        for row in array:
            for value in row:
                if value and abs(value) > eps:
                    n += 1
        return n
    
    @staticmethod
    def FillHistogram( array, title = 'array', xlabel = 'x (mm)', ylabel = 'y (mm)', zlabel = 'Q (pes)' ):
        '''
            Fill an histogram with the data of the array.
        '''
        hxy = NEXT.TrackingPlane.GetHistogram( title, title )
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
        return numpy.array( [ [ histogram.GetBinContent( i, j ) for j in range( 1, TrackingPlane.Nbins + 1 ) ]
                                                                for i in range( 1, TrackingPlane.Nbins + 1 ) ] )



