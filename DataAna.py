from FourierTransforms import Transformer,Shrink
from ROOT import TCanvas, TPad, gStyle
from shelve import open as openshelve
from Plots import MakeGif
from TrackBuilder import *
from Tools import Arrays


gStyle.SetOptStat('')

datadir  = '/Users/Gonzalo/Desktop/Tracks/'
datafile = datadir + 'CsRTracks.shelve'
datafile = datadir + 'TlRTracks.shelve'
output = './dump/'

Nevts = int(100)
data = openshelve( datafile )

for i in map( str, range(1000) ):
    print i
    track = data.get( i )
    if track.E < 20000 or track.E > 30000:
        continue
    print track
    raw_input()
#    if track.E <  11160. or track.E > 11640. or track.nslices < 40:
#        continue
    isources = []
    fsources = []
    signals  = []
    
    for j,slice in enumerate(track.slices):
        if slice.nsipms < 4:
            continue
        slice.ComputeSiPMsMatrix()
        algorithm = Transformer( Shrink(slice.SiPMsMatrix) )
        algorithm.Compute(cut = 0.94)
        
        isources.append( algorithm.sources[0] )
        fsources.append( algorithm.fsource    )
        signals .append( slice.SiPMsMatrix    )

    nslices = len(isources)
    if not nslices:
        continue

    imax = max( map( lambda x: Arrays.Max(x)[2], isources ) )
    fmax = max( map( lambda x: Arrays.Max(x)[2], fsources ) )


    for j in range(nslices):
        isources[j] = Arrays.FillHistogram( isources[j], 'Fourier method'  , name = 'FM' + i + str(j) )
        fsources[j] = Arrays.FillHistogram( fsources[j], 'Iterative method', name = 'IM' + i + str(j) )
        signals [j] = Arrays.FillHistogram( signals [j], 'SiPMs signal'    , name = 'S'  + i + str(j) ).RebinX(5).RebinY(5)
        isources[j].SetMaximum(imax)
        fsources[j].SetMaximum(fmax)

        c  = TCanvas()
        p1 = TPad( 'a', 'a', 0.00, 0.05, 0.66, 0.95 ); p1.SetGrid(); p1.Draw()
        p2 = TPad( 'b', 'b', 0.67, 0.50, 1.00, 1.00 ); p2.SetGrid(); p2.Draw()
        p3 = TPad( 'c', 'c', 0.67, 0.00, 1.00, 0.50 ); p3.SetGrid(); p3.Draw()
        p1.cd(); fsources[j].Draw('zcol')
        p2.cd(); isources[j].Draw('zcol')
        p3.cd(); signals [j].Draw('zcol')

        c.SaveAs( 'dump/track{0}slice{1}.png'.format(i,j) )

    MakeGif( ['track{0}slice{1}'.format(i,j) for j in range(nslices) ], 'dump/', output = './dump/track{0}'.format(i) )
