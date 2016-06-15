'''
    Analysis of the Fourier toy MC.
    
    @Author: G. Martinez Lema
    @Date  : 22/10/2014
'''

from __future__ import division
from FourierToy import *

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

def ChargeStudies( Nexp = 1e3 ):
    R  = 1/50
    N0 = 1e5
    N  = N0 * R
    hratio  = ROOT.TH1D( 'Q/QT'     , '', 100, 0, R )
    hQvsQT  = ROOT.TH2D( 'QvsQT'    , '', 100, 0, N0, 100, 0, N )
    hQvscls = ROOT.TH2D( 'Qvscls'   , '',   7, 0, 7, 100, 0, N )
    hrvscls = ROOT.TH2D( 'ratio/QT' , '',   7, 0, 7, 100, 0, R )
    
    random = ROOT.TRandom3(0)
    Nexp = int( Nexp )
    for i in range(Nexp):
        Nph = random.Uniform( 1e3, N0 )
        x0  = random.Uniform(  -5, 5 )
        y0  = random.Uniform(  -5, 5 )
        cls = math.sqrt( ( abs(x0) - 1 )**2 + ( abs(y0) - 1 )**2 )
        
        simulation = Simulator( Pointlike(x0,y0), Nph, 5. )
        xy_full = Tools.Arrays.FromHistogram( simulation.GetXYDistribution( FullInformation = True ) )
        xy_disc = NEXT.TrackingPlane.Discretize( xy_full )
        
        Nfull = Tools.Arrays.Sum( xy_full )
        Ndisc = Tools.Arrays.Sum( xy_disc )
        Ratio = Ndisc/Nfull

        hratio .Fill( Ratio )
        hQvsQT .Fill( Nfull , Ndisc )
        hQvscls.Fill(    cls, Ndisc )
        hrvscls.Fill(    cls, Ratio )
        del simulation

    return hratio, hQvsQT, hQvscls, hrvscls

import Plots
import time
print 'running'
t0 = time.time()

hratio, hQvsQT, hQvscls, hrvscls = ChargeStudies( 10000 )
hratio.GetXaxis().SetTitle('Q_{SiPMs}/Q_{Total}');hratio.GetYaxis().SetTitle('Entries')
hQvsQT.GetXaxis().SetTitle('Q_{Total}');hQvsQT.GetYaxis().SetTitle('Q_{SiPMs}')
hQvscls.GetXaxis().SetTitle('Distance to closest SiPM (mm)');hQvscls.GetYaxis().SetTitle('Q_{SiPMs}')
hrvscls.GetXaxis().SetTitle('Distance to closest SiPM (mm)');hrvscls.GetYaxis().SetTitle('Q_{SiPMs}/Q_{Total}')

histos = [ hratio, hrvscls, hQvsQT, hQvscls.ProfileX() ]
c = Plots.PutInCanvas( histos, ['','','zcol',''] )
c.cd(2)
histos[1].ProfileX().Draw('same')

#sim = Simulator( UniformCircle( 4., 5., 5.), 1e7, 0. )
#xy  = sim.GetXYDistribution( FullInformation = True )
#xy  = Tools.Arrays.FromHistogram( xy )
#xy  = NEXT.TrackingPlane.Discretize( xy )
#psf = Tools.Arrays.MakePSF( xy )
#psf.SetMarkerStyle(20)
#psf.SetMarkerSize(1)
#fun = ROOT.TF1('fit','[0]/( 1 + [1]*x^2 )^1.5')
#fun.SetParameters( 1,1e-3 )
#fun.SetParLimits(0,0,10)
#fun.SetParLimits(1,0,10)
#psf.Fit( fun )
#psf.Draw('AP')

print time.time() - t0