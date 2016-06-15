from ROOT import *

file = TFile( 'str.root' )

#plot = TH2D( file.Get( 'Pointlike_plot;1' ) )
#xy   = TH2D( file.Get( 'Pointlike_xy;1'   ) )
#r    = TH1D( file.Get( 'Pointlike_r;1'    ) )

plotplus = [ TH2D( file.Get( 'Pointlike_plot;' + str(i) ) ) for i in range(1,6) ]
xyplus   = [ TH2D( file.Get( 'Pointlike_xy;'   + str(i) ) ) for i in range(1,6) ]
rplus    = [ TH1D( file.Get( 'Pointlike_r;'    + str(i) ) ) for i in range(1,6) ]

plotplus = [ TH2D( file.Get( 'Pointlike_plot;' + str(i) ) ) for i in range(1,6) ]
xyplus   = [ TH2D( file.Get( 'Pointlike_xy;'   + str(i) ) ) for i in range(1,6) ]
rplus    = [ TH1D( file.Get( 'Pointlike_r;'    + str(i) ) ) for i in range(1,6) ]

plot = plotplus[0]
xy   =   xyplus[0]
r    =    rplus[0]

for i in range(1,5):
    plotplus[0].Add( plotplus[i] )
    xyplus[0]  .Add(   xyplus[i] )
    rplus[0]   .Add(    rplus[i] )

fileout = TFile( 'str_out.root', 'recreate' )
plot.Write()
xy  .Write()
r   .Write()

