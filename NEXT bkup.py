'''
    Module with information about NEXT.
    
    Author: G. Martinez Lema
    Date  : 17/10/2014
    
'''

import ROOT
import numpy
import math
import copy

class TrackingPlane:
    '''
        Tracking plane data.
    '''
    Pitch   = 1.
    Range   = 75.
    Nbins   = int( 2 * Range ) + 1
    NbinsS  = Nbins - 1
    NSiPMs  = 256
    NSiPMs_side = int( round(  math.sqrt(NSiPMs) ) )
    
    Bins    = [ -Range + 1.*i for i in range( Nbins ) ]
    BinsS   = Bins[:-1]
    SiPMs   = [ -Range + 10*i for i in range( NSiPMs_side ) ]
    SiPMs_indices = [ int( x + Range ) for x in SiPMs ]
    SiPMsS_indices = SiPMs_indices[:-1]

    @staticmethod
    def GetHistogram( name, title = '' ):
        '''
            Return an empty histogram with the granularity of the tracking plane.
        '''
        
        if not title:
            title = name
        return ROOT.TH2D( name, title,
                          TrackingPlane.Nbins, -(TrackingPlane.Range + .5), (TrackingPlane.Range + .5),
                          TrackingPlane.Nbins, -(TrackingPlane.Range + .5), (TrackingPlane.Range + .5))
    
    @staticmethod
    def Discretize( arg ):
        '''
            Erase everything but the information of the SiPMs from the argument.
        '''
        
        if isinstance( arg, ROOT.TH2 ):
            return TrackingPlane._DiscretizeHistogram( arg )
        else:
            return TrackingPlane._DiscretizeArray( arg )
    
    @staticmethod
    def _DiscretizeHistogram( hxy ):
        '''
            Erase everything but the information of the SiPMs from the histogram.
        '''
        
        hxy_SiPMs = copy.copy( hxy )
        for i in range( 1, TrackingPlane.Nbins + 1 ):
            for j in range( 1, TrackingPlane.Nbins + 1 ):
                if (hxy_SiPMs.GetXaxis().GetBinCenter(i) in TrackingPlane.SiPMs and
                    hxy_SiPMs.GetYaxis().GetBinCenter(j) in TrackingPlane.SiPMs):
                    continue
                bin  = hxy_SiPMs.GetBin(   i,  j  )
                hxy_SiPMs.SetBinContent( bin,  0  )
                hxy_SiPMs.SetBinError  ( bin,  0  )
        
        return hxy_SiPMs
    
    @staticmethod
    def _DiscretizeArray( array ):
        '''
            Erase everything but the information of the SiPMs from the array.
        '''
        
        Darray = copy.deepcopy( array )
        rows, cols = array.shape
        for i in range( rows ):
            for j in range( cols ):
                if (i - TrackingPlane.Range in TrackingPlane.SiPMs and
                    j - TrackingPlane.Range in TrackingPlane.SiPMs):
                    continue
                Darray[i][j] = 0.
        return Darray

    
__all__ = [ 'TrackingPlane' ]
