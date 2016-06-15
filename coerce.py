from ROOT import *
from DatabaseNEW import *
from array import array

R = TRandom3(0)
r = 222.5
r2 = r**2

'''
for i in xrange(int(1e7)):
    x = R.Uniform(-r,r)
    y = R.Uniform(-r,r)
    if (x**2+y**2) > r2: continue
    
''' 

#h = TH2F('a','',441,-220.5,220.5,441,-220.5,220.5)
h = TH2F('a','',445,-222.5,222.5,445,-222.5,222.5)
g = TGraph()
for i in xrange(int(1e6)):
    x = R.Gaus(0,1)
    y = R.Gaus(0,1)
    n = (x**2+y**2)**0.5
    g.SetPoint( i, x/n*r, y/n*r )

for id,xy in SiPM_map.items(): h.SetBinContent(h.Fill(*xy),id)

g.Draw('ap')
h.Draw('zcolsame')
raw_input()
