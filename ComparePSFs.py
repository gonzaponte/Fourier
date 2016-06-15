from Statistics import *
from Plots import *

PSF = lambda p0,p1,p2,p3: lambda x: p0/1e3 * ( 1 + p1/1e3 * x**2 + p2/1e6 * x**4 )**-1.5 + p3/1e3


F_RXall          = PSF( 42.78, 1.089, 1.123, 0.8029 )
F_RXcenter       = PSF( 70.58, 1.999, 4.181, 1.5920 )
F_RXoverSiPM     = PSF( 66.28, 2.241, 2.066, 0.9403 )
F_BlobBunched    = PSF( 52.98, 2.542, 1.101, 1.6590 )
F_BlobIndividual = PSF( 53.79, 2.535, 1.199, 1.6910 )
F_Mip            = PSF( 71.69, 2.380, 2.406, 1.0000 )

F_G2             = PSF( 496.2, 29.43, 21.42, 0.7888 )
F_G4             = PSF( 320.0, 12.40, 26.52, 0.3858 )
F_G6             = PSF( 202.5, 6.086, 14.58, 0.2722 )
F_G8             = PSF( 139.6, 3.645, 7.932, 0.2289 )
F_G10            = PSF( 99.08, 2.381, 4.047, 0.1339 )

F_C2_8           = PSF( 134.8, 3.346, 7.842, 0.2321 )
F_C2_10          = PSF( 97.15, 2.262, 4.025, 0.1159 )
F_C5_8           = PSF( 130.6, 3.499, 6.444, 0.1801 )
F_C5_10          = PSF( 94.42, 2.242, 3.628, 0.1132 )

for b in [True,False]:
    RXall          = Distribution( F_RXall         , 0, 50, 100, normalized = b )
    RXcenter       = Distribution( F_RXcenter      , 0, 50, 100, normalized = b )
    RXoverSiPM     = Distribution( F_RXcenter      , 0, 50, 100, normalized = b )
    BlobBunched    = Distribution( F_BlobBunched   , 0, 50, 100, normalized = b )
    BlobIndividual = Distribution( F_BlobIndividual, 0, 50, 100, normalized = b )
    Mip            = Distribution( F_Mip           , 0, 50, 100, normalized = b )

    G2             = Distribution( F_G2            , 0, 50, 100, normalized = b )
    G4             = Distribution( F_G4            , 0, 50, 100, normalized = b )
    G6             = Distribution( F_G6            , 0, 50, 100, normalized = b )
    G8             = Distribution( F_G8            , 0, 50, 100, normalized = b )
    G10            = Distribution( F_G10           , 0, 50, 100, normalized = b )

    C2_8           = Distribution( F_C2_8          , 0, 50, 100, normalized = b )
    C2_10          = Distribution( F_C2_10         , 0, 50, 100, normalized = b )
    C5_8           = Distribution( F_C5_8          , 0, 50, 100, normalized = b )
    C5_10          = Distribution( F_C5_10         , 0, 50, 100, normalized = b )


    print 'RXall          integral = ', RXall         .Integral()
    print 'RXcenter       integral = ', RXcenter      .Integral()
    print 'RXoverSiPM     integral = ', RXoverSiPM    .Integral()
    print 'BlobBunched    integral = ', BlobBunched   .Integral()
    print 'BlobIndividual integral = ', BlobIndividual.Integral()
    print 'Mip            integral = ', Mip           .Integral()
    print ''
    print 'Gauss 2        integral = ', G2            .Integral()
    print 'Gauss 4        integral = ', G4            .Integral()
    print 'Gauss 6        integral = ', G6            .Integral()
    print 'Gauss 8        integral = ', G8            .Integral()
    print 'Gauss 10       integral = ', G10           .Integral()
    print ''
    print 'Circle 2 + 8   integral = ', C2_8          .Integral()
    print 'Circle 2 + 10  integral = ', C2_10         .Integral()
    print 'Circle 5 + 8   integral = ', C5_8          .Integral()
    print 'Circle 5 + 10  integral = ', C5_10         .Integral()
    print ''


HRXall          = RXall         .PDFHistogram()
HRXcenter       = RXcenter      .PDFHistogram()
HRXoverSiPM     = RXoverSiPM    .PDFHistogram()
HBlobBunched    = BlobBunched   .PDFHistogram()
HBlobIndividual = BlobIndividual.PDFHistogram()
HMip            = Mip           .PDFHistogram()

HG2             = G2            .PDFHistogram()
HG4             = G4            .PDFHistogram()
HG6             = G6            .PDFHistogram()
HG8             = G8            .PDFHistogram()
HG10            = G10           .PDFHistogram()

HC2_8           = C2_8          .PDFHistogram()
HC2_10          = C2_10         .PDFHistogram()
HC5_8           = C5_8          .PDFHistogram()
HC5_10          = C5_10         .PDFHistogram()

HRXall         .SetLineWidth(3)
HRXcenter      .SetLineWidth(3)
HRXoverSiPM    .SetLineWidth(3)
HBlobBunched   .SetLineWidth(3)
HBlobIndividual.SetLineWidth(3)
HMip           .SetLineWidth(3)

HG2            .SetLineWidth(3)
HG4            .SetLineWidth(3)
HG6            .SetLineWidth(3)
HG8            .SetLineWidth(3)
HG10           .SetLineWidth(3)

HC2_8          .SetLineWidth(3)
HC2_10         .SetLineWidth(3)
HC5_8          .SetLineWidth(3)
HC5_10         .SetLineWidth(3)

#HRXall         .SetLineStyle(1)
#HRXcenter      .SetLineStyle(2)
#HRXoverSiPM    .SetLineStyle(6)
#HBlobBunched   .SetLineStyle(8)
#HBlobIndividual.SetLineStyle(9)
HMip           .SetLineStyle(4)
#
#HG10           .SetLineStyle(1)
#HG10           .SetLineStyle(2)
#HG10           .SetLineStyle(3)
#HG10           .SetLineStyle(4)
#HG10           .SetLineStyle(5)
#
#HC2_8          .SetLineStyle(6)
#HC2_10         .SetLineStyle(7)
#HC5_8          .SetLineStyle(8)
#HC5_10         .SetLineStyle(9)


ROOT.gStyle.SetOptStat('')

c = Merge1D([
#              HRXall         ,
#              HRXcenter      ,
#              HRXoverSiPM    ,
#              HBlobBunched   ,
#              HBlobIndividual,
              HMip           ,
              HG2            ,
              HG4            ,
              HG6            ,
              HG8            ,
              HG10           ,
             
              HC2_8          ,
              HC2_10         ,
              HC5_8          ,
              HC5_10         ],
              '#Deltar (mm)', 'Q_{SiPM} / Q_{Total}' )


legend = ROOT.TLegend(0.6,0.6,0.9,0.9)

#legend.AddEntry( HRXall         , 'RX all'          , 'L' )
#legend.AddEntry( HRXcenter      , 'RX center'       , 'L' )
#legend.AddEntry( HRXoverSiPM    , 'RX over SiPM'    , 'L' )
#legend.AddEntry( HBlobBunched   , 'Blob bunched'    , 'L' )
#legend.AddEntry( HBlobIndividual, 'Blob individual' , 'L' )
legend.AddEntry( HMip           , 'Mip'             , 'L' )

legend.AddEntry( HG2            , 'Gauss 2  mm'     , 'L' )
legend.AddEntry( HG4            , 'Gauss 4  mm'     , 'L' )
legend.AddEntry( HG6            , 'Gauss 6  mm'     , 'L' )
legend.AddEntry( HG8            , 'Gauss 8  mm'     , 'L' )
legend.AddEntry( HG10           , 'Gauss 10 mm'     , 'L' )

legend.AddEntry( HC2_8          , 'Circle 2 + 8 mm' , 'L' )
legend.AddEntry( HC2_10         , 'Circle 2 + 10 mm', 'L' )
legend.AddEntry( HC5_8          , 'Circle 5 + 8 mm' , 'L' )
legend.AddEntry( HC5_10         , 'Circle 5 + 10 mm', 'L' )

legend.Draw('same')


#[0]/(1+exp((x-[1])/[2]))+[3]
#Fit to #frac{p0}{1+exp(#frac{x-p1}{p2})} + p3

