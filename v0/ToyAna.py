'''
    Analysis of the Fourier toy MC.
    
    @Author: G. Martinez Lema
    @Date  : 22/10/2014
'''


from FourierToy import *
import time

t0 = time.time()




xy_diffusion = 0. # mm
outputdir = './Output/'
outputfile = outputdir + 'structures_{0}mm.root'.format(int(xy_diffusion))

def ComputePointlike( N = 1e8, fname = outputfile ):
    print 'Computing pointlike source with {0:<1.1e} photons'.format( N )
    pointlike = Simulator( Pointlike(), N, xy_diffusion )
    pointlike.SaveToROOTFile( fname, 'update' )

def ComputeUniformCircle( N = 1e8, fname = outputfile ):
    print 'Computing uniform circle source with {0:<1.1e} photons'.format( N )
    ucircle = Simulator( UniformCircle(4.), N, xy_diffusion )
    ucircle.SaveToROOTFile( fname, 'update' )

def ComputeStraightLine( N = 1e8, fname = outputfile ):
    print 'Computing straight line source with {0:<1.1e} photons'.format( N )
    strline = Simulator( StraightLine( -5, -5, 10 * 2**.5, math.pi/4), N, xy_diffusion )
    strline.SaveToROOTFile( fname, 'update' )

def MakePSF( struct, fname, diffusion = xy_diffusion ):
    print fname
    file = ROOT.TFile( fname, 'recreate' )
    
    simu = Simulator( struct, 1e7, diffusion )
    xy  = simu.GetXYDistribution()
    xy  = Tools.Arrays.FromHistogram( xy )
    PSF = Tools.Arrays.MakePSF( xy )
    PSF.SetMarkerStyle(20)
    PSF.SetMarkerSize(1)
    
    psf_perfect = ROOT.TF1('psf_perfect', '[0]/( 1 + [1]*x^2           )^1.5'       )
    psf_bg      = ROOT.TF1('psf_bg'     , '[0]/( 1 + [1]*x^2           )^1.5 + [2]' )
    psf_r       = ROOT.TF1('psf_r'      , '[0]/( 1 + [1]*x^2 + [2]*x^4 )^1.5'       )
    psf_bgr     = ROOT.TF1('psf_bgr'    , '[0]/( 1 + [1]*x^2 + [2]*x^4 )^1.5 + [3]' )
    psf_perfect.SetParameters( 1,1e-3             )
    psf_bg     .SetParameters( 1,1e-3, 1e-6       )
    psf_r      .SetParameters( 1,1e-3, 1e-6       )
    psf_bgr    .SetParameters( 1,1e-3, 1e-6, 1e-6 )
    psf_perfect.SetParLimits(0,0,10); psf_perfect.SetParLimits(1,0,1)
    psf_bg     .SetParLimits(0,0,10); psf_bg     .SetParLimits(1,0,1); psf_bg .SetParLimits(2,0,1)
    psf_r      .SetParLimits(0,0,10); psf_r      .SetParLimits(1,0,1); psf_r  .SetParLimits(2,0,1)
    psf_bgr    .SetParLimits(0,0,10); psf_bgr    .SetParLimits(1,0,1); psf_bgr.SetParLimits(2,0,1); psf_bgr.SetParLimits(3,0,1)

    PSF        .Write()
    psf_perfect.Write()
    psf_bg     .Write()
    psf_r      .Write()
    psf_bgr    .Write()

    file.Close()


N = 1e9
#ComputePointlike(N)
#ComputeUniformCircle(N)
#ComputeStraightLine(N)
#Q_{relative}
#Deltar (mm)
#[0]/(1+exp((x-[1])/[2]))
#Fit to #frac{p0}{1+exp(#frac{x-p1}{p2})}

#MakePSF( Pointlike(4.,4.), 'pointlikepsf.root', 0. )
#MakePSF( UniformCircle(20.,4.,4.), 'circlepsf_20mm.root', 0. )
#MakePSF( Pointlike(4.,4.), 'gaussianpsf_5mm.root', 5. )
#MakePSF( Pointlike(4.,4.), 'gaussianpsf_10mm.root', 10. )
#MakePSF( Pointlike(4.,4.), 'gaussianpsf_20mm.root', 20. )
MakePSF( Pointlike(4.,4.), 'gaussianpsf_2mm.root', 2. )
MakePSF( Pointlike(4.,4.), 'gaussianpsf_4mm.root', 4. )
MakePSF( Pointlike(4.,4.), 'gaussianpsf_6mm.root', 6. )
MakePSF( Pointlike(4.,4.), 'gaussianpsf_8mm.root', 8. )
#MakePSF( UniformCircle(5.,4.,4.), 'circlepsf_5mm_8mm.root', 8. )
#MakePSF( UniformCircle(5.,4.,4.), 'circlepsf_5mm_10mm.root', 10. )
#MakePSF( UniformCircle(2.,4.,4.), 'circlepsf_2mm_8mm.root', 8. )
#MakePSF( UniformCircle(2.,4.,4.), 'circlepsf_2mm_10mm.root', 10. )

dt  = time.time() - t0
h   = dt // 3600
dt -= h * 3600
min = dt // 60
s   = dt - min * 60

print '{0} h {1} min {2} s'.format( h, min, s )



