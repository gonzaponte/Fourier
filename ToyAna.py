'''
    Analysis of the Fourier toy MC.
    
    @Author: G. Martinez Lema
    @Date  : 22/10/2014
'''


from FourierToy import *

xy_diffusion = 4. # mm
outputdir = './Output/'
outputfile = outputdir + 'structures.root'

def ComputePointlike( N = 1e8, fname = outputfile ):
    pointlike = Simulator( Pointlike(), N )
    pointlike.SaveToROOTFile( fname, 'update' )

def ComputeUniformCircle( N = 1e8, fname = outputfile ):
    ucircle = Simulator( UniformCircle(4.), N )
    ucircle.SaveToROOTFile( fname, 'update' )

def ComputeStraightLine( N = 1e8, fname = outputfile ):
    strline = Simulator( StraightLine( -5, -5, 10 * 2**.5, math.pi/4), N )
    strline.SaveToROOTFile( fname, 'update' )


N = 1e8
ComputePointlike(N)
ComputeUniformCircle(N)
ComputeStraightLine(N)




