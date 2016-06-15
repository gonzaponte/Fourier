from DatabaseNEW import SiPM_map
from ROOT import TH2F,TCanvas,gStyle

gStyle.SetOptStat('')

hboard = TH2F( 'NEWboards', 'NEW SiPMs ID map: board number', 500, -250, 250, 500, -250, 250 )
hid    = TH2F( 'NEWSiPMs' , 'NEW SiPMs ID map: SiPM number',    8,    0,  80,   8,    0,  80 )

for i in range(30):
    try:
        x,y = SiPM_map[ i*1000 + 28 ]
        for j in range(i):
            hboard.Fill( x, y )
    except:
        pass

for i in range(64):
    x,y = SiPM_map[ 11000 + i ]
    [ hid.Fill( x, y ) for j in range(i) ]

hboard.Scale(1000)

c = TCanvas()
c.Divide(2,1)
c.cd(1)
hboard.Draw('text')
c.cd(2)
hid.Draw('text0')