'''
    A toy MC to simulate the anode signal produced by different structures that emmit light isotropically from the EL mesh.
    
    @Author: G. Martinez Lema
    @Date  : 16/10/2014
'''

import ROOT
import math
import copy
import NEXT
import Tools

class Structure:
    '''
        The 2D/3D-shape that produces light at the EL mesh.
    '''
    
    def __init__( self, name ):
        '''
            Constructor with the name of the structure.
        '''
        self.name   = name
        self.z0     = 4.5
        self.random = ROOT.TRandom3(0)

    def GetPointXYZ( self ):
        '''
            Returns a random point within the structure. Must be implemented in each particular case.
        '''
        return 0., 0., 4.5

    def Plot( self, Npoints = 1e6, sigma_diffusion = 0. ):
        '''
            Build a 2D-histogram with the shape of the structure at the EL mesh (z-integrated) with the tracking plane granularity.
        '''
        Npoints = int(Npoints)
        hxy_structure = NEXT.TrackingPlane.GetHistogram( self.name + '_plot', self.name )
        
        for i in xrange( Npoints ):
            x, y, z = self.GetPointXYZ()
            
            if sigma_diffusion: #diffusion
                x += self.random.Gaus( 0, sigma_diffusion )
                y += self.random.Gaus( 0, sigma_diffusion )
            
            hxy_structure.Fill( x, y )

        return hxy_structure


class Pointlike( Structure ):
    '''
        Pointlike source.
    '''
    def __init__( self, x0 = 0., y0 = 0. ):
        Structure.__init__( self, 'Pointlike' )
        self.x0 = x0
        self.y0 = y0

    def GetPointXYZ( self ):
        return self.x0, self.y0, self.z0

class UniformCircle( Structure ):
    '''
        Uniform circle.
    '''
    def __init__( self, r = 1., x0 = 0., y0 = 0.):
        Structure.__init__( self, 'Uniform circle' )
        self.x0 = x0
        self.y0 = y0
        self.r  = r

    def GetPointXYZ( self ):
        while True:
            x = self.random.Uniform( -1, 1 )
            y = self.random.Uniform( -1, 1 )
            if x**2 + y**2 <= 1:
                return x * self.r + self.x0, y * self.r + self.y0, self.z0

class GaussianCircle( Structure ):
    '''
        Gaussian circle.
    '''
    def __init__( self, sigma = 1., x0 = 0., y0 = 0., ):
        Structure.__init__( self, 'Gaussian circle' )
        self.x0 = x0
        self.y0 = y0
        self.sigma = sigma
    
    def GetPointXYZ( self ):
        x   = self.random.Gaus( 0., self.sigma  )
        y   = self.random.Gaus( 0., self.sigma  )
        return x + self.x0, y + self.y0, self.z0

class UniformEllipse( Structure ):
    '''
        Uniform ellipse.
    '''
    def __init__( self, a = 2., b = 1., x0 = 0., y0 = 0. ):
        Structure.__init__( self, 'Ellipse' )
        self.a = a
        self.b = b
        self.x0 = x0
        self.y0 = y0

    def GetPointXYZ( self ):
        while True:
            x = self.random.Uniform( -1, 1 )
            y = self.random.Uniform( -1, 1 )
            if x**2 + y**2 <= 1:
                return x * self.a + self.x0, y * self.b + self.y0, self.z0

class GaussianEllipse( Structure ):
    '''
        Gaussian ellipse.
    '''
    def __init__( self, sigmax = 3., sigmay = 1., x0 = 0., y0 = 0., ):
        Structure.__init__( self, 'Gaussian ellipse' )
        self.x0 = x0
        self.y0 = y0
        self.sigmax = sigmax
        self.sigmay = sigmay
    
    def GetPointXYZ( self ):
        x   = self.random.Gaus( 0., self.sigmax  )
        y   = self.random.Gaus( 0., self.sigmay  )
        return x + self.x0, y + self.y0, self.z0


class StraightLine( Structure ):
    '''
        A straight line in the x-y plane.
    '''
    def __init__( self, x0 = 0., y0 = 0., L = 10., theta = math.pi/4. ):
        Structure.__init__( self, 'Straight line' )
        self.x0 = x0
        self.y0 = y0
        self.L  = L
        self.theta = theta
        
        self.cos = math.cos( self.theta )
        self.sin = math.sin( self.theta )

    def GetPointXYZ( self ):
        r = self.random.Uniform( 0., self.L )
        return r * self.cos + self.x0, r * self.sin + self.y0, self.z0

class Simulator:
    '''
        The master class that produces results.
    '''
    
    def __init__( self, structure, Nphotons = 1e8, diffusion = 0. ):
        '''
            Constructor. Give it some structure.
        '''
        self.structure = structure
        self.random    = ROOT.TRandom3(0)
        self.Nphotons  = int(Nphotons)
        self.diffusion = diffusion
        self.hxy       = NEXT.TrackingPlane.GetHistogram( self.structure.name + '_xy', self.structure.name )
        self.hr        = ROOT.TH1D( self.structure.name + '_r' , self.structure.name,
                                    1000, 0, NEXT.TrackingPlane.Range * math.sqrt(2) )

        self._Generate()

    def _Generate( self ):
        '''
            Generate the signal in the SiPMs.
        '''
        for i in xrange( self.Nphotons ):
            x0, y0, z0  = self.structure.GetPointXYZ()
            
            if self.diffusion:
                x0 += self.random.Gaus( 0, self.diffusion )
                y0 += self.random.Gaus( 0, self.diffusion )
            
            theta = math.acos( self.random.Uniform() )
            phi   = 2 * math.pi * self.random.Uniform()

            x  = z0 * math.tan( theta ) * math.cos( phi )
            y  = x  * math.tan( phi )
            x += x0
            y += y0
            r  = math.sqrt( x**2 + y**2 )

            self.hxy.Fill( x, y )
            self.hr .Fill( r )

    def GetXYDistribution( self, FullInformation = False ):
        '''
            Return the 2D histogram. Set FullInformation to true to get the signal in the whole plane.
        '''
        return self.hxy if FullInformation else NEXT.TrackingPlane.Discretize( self.hxy )
    
    def GetRDistribution( self ):
        '''
            Return the 1D histogram.
        '''
        return self.hr
    
    def GetStructure( self, Npoints = 1e7 ):
        '''
            Return the histogram of the structure.
        '''
        return self.structure.Plot( Npoints, self.diffusion )

    def SaveToROOTFile( self, rootfile, option = 'recreate' ):
        '''
            Save histograms to a given ROOT file. The input can be either the file or its name. In the latter case, open it with the specified option.
        '''
        if isinstance( rootfile, str ):
            rootfile = ROOT.TFile( rootfile, option )
        
        self.GetStructure().Write()
        self.hxy.Write()
        self.hr .Write()
        rootfile.Close()

if __name__ == '__main__':
    ### Some examples with few points to perform rapidly
    ### Comment everything but the structure to be visualized
    
    structure = Pointlike( x0 = 15., y0 = 40. )
    structure = UniformCircle( r = 5. )
    structure = GaussianCircle( sigma = 25. )
    structure = GaussianEllipse( sigmax = 30, sigmay = 10 )
    structure = StraightLine( L = 20. )

    simu = Simulator( structure, Nphotons = 1e5 )
    hxy  = simu.GetXYDistribution( FullInformation = True )
    hstr = simu.GetStructure( Npoints = 1e5 )

    c = ROOT.TCanvas()
    c.Divide(1,2)
    c.cd(1)
    hstr.Draw('zcol')
    c.cd(2)
    hxy .Draw('zcol')




