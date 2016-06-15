'''
  TrackReader.py

  Description:
  Reads the trees data and generates the classes.

  Author:  Gonzalo
  Date:   19/01/2014

  Last update 05/02/2014 <- Gonzalo
'''



import ROOT
from math import *
from DatabaseMC import *
from Utils import *
from reader import IReader

#Track instance
class track:
    def __init__(self,E,Q,nslices,Emaxpos,Qmaxpos,isright,slices):
        self.E       = E
        self.Q       = Q
        self.nslices = nslices
        self.Emaxpos = Emaxpos
        self.Qmaxpos = Qmaxpos
        self.isright = isright
        self.slices  = slices
    
    def Emax(self):
        ''' Energy of the most energetic slice.'''
        self.Emax = self.slices[self.Emaxpos].E

    def EQmax(self):
        ''' Energy of the most charged slice.'''
        self.EQmax = self.slices[self.Qmaxpos].E

    def Qmax(self):
        ''' Charge of the most charged slice.'''
        self.Qmax = self.slices[self.Qmaxpos].Q

    def QEmax(self):
        ''' Charge of the most energetic slice.'''
        self.QEmax = self.slices[self.Emaxpos].Q

    def Xmax(self):
        ''' X of the most energetic slice.'''
        self.Xmax = self.slices[self.Emaxpos].X

    def Ymax(self):
        ''' Y of the most energetic slice.'''
        self.Ymax = self.slices[self.Emaxpos].Y

    def Zmax(self):
        ''' Z of the most energetic slice.'''
        self.Zmax = self.slices[self.Emaxpos].Z

    def GetAllE( self ):
        ''' Return a list with the energy of each slice.'''
        return [ slice.E for slice in self.slices ]

    def GetAllEcorrected( self ):
        ''' Return a list with the corrected energy of each slice.'''
        return [ slice.Ecor for slice in self.slices ]

    def GetAllQ( self ):
        ''' Return a list with the charge of each slice.'''
        return [ slice.Q for slice in self.slices ]

    def GetAllX( self ):
        ''' Return a list with the x position of each slice.'''
        return [ slice.X for slice in self.slices ]

    def GetAllY( self ):
        ''' Return a list with the y position of each slice.'''
        return [ slice.Y for slice in self.slices ]

    def GetAllZ( self ):
        ''' Return a list with the z position of each slice.'''
        return [ slice.Z for slice in self.slices ]

    def GetAllIDmax( self ):
        ''' Return a list with the ID of the SiPM with greatest signal of each slice.'''
        return [ slice.E for slice in self.slices ]

    def GetAllNsipms( self ):
        ''' Return a list with the number of sipms touched of each slice.'''
        return [ slice.nsipms for slice in self.slices ]

    def GetAllSiPMs( self ):
        ''' Return a list with the SiPMs of each slice.'''
        return [ slice.sipms for slice in self.slices ]

    def GetAllPMTs( self ):
        ''' Return a list with the PMTs of each slice.'''
        return [ slice.pmts for slice in self.slices ]

    def ComputeEcorrected( self ):
        ''' Return the energy of the track with applied corrections.'''
        self.Ecor = sum( self.GetAllEcorrected() )
    
    def __repr__(self):
        return 'E: {0}   Q: {1}   nslices: {2}\nSlices:\n{3}'.format( self.E, self.Q, self.nslices, '\n'.join( map( str, self.slices ) ) )


#Slice instance
class slice:
    def __init__(self,E,Q,X,Y,Z,IDmax,nsipms,sipms,pmts):
        self.E      = E
        self.Ecor   = E
        self.Q      = Q
        self.Qcor   = Q
        self.X      = X
        self.Y      = Y
        self.R      = sqrt( X**2 + Y**2 )
        self.phi    = Angle( X, Y )
        self.corona = Corona( self.R, R0corona )
        self.sector = PhiSector( self.phi, Phi0corona )
        self.Z      = Z
        self.IDmax  = IDmax
        self.nsipms = nsipms
        self.sipms  = sipms
        self.pmts   = pmts
    
    def Qmax(self):
        ''' Return the charge of the sipm with greatest signal.'''
        self.Qmax = self.sipms[0].Q if self.sipms else 0.
    
    def Xmax(self):
        ''' Return the x position of the sipm with greatest signal.'''
        self.Xmax = SiPM_map[self.IDmax][0]

    def Ymax(self):
        ''' Return the y position of the sipm with greatest signal.'''
        self.Ymax = SiPM_map[self.IDmax][1]

    def Rmax(self):
        ''' Return the r position of the sipm with greatest signal.'''
        self.Rmax = sqrt(self.Xmax()**2 + self.Ymax()**2)
    
    def GetSiPM(self,ID):
        ''' Return the SiPM with the specified ID.'''
        for sipm in self.sipms:
            if sipm.ID == ID:
                return sipm
        print 'ID {0} not found in slice.GetSiPM'.format(ID)
        return False
        
    def GetPMT(self,ID):
        ''' Return the PMT with the specified ID.'''
        return self.pmts[ID]

    def GetAllE( self ):
        ''' Return a list with the signal of each PMT.'''
        return [ pmt.E for pmt in self.pmts.values() ]

    def GetAllQ( self ):
        ''' Return a list with the signal of each SiPM.'''
        return [ sipm.Q for sipm in self.sipms ]

    def GetAllX( self ):
        ''' Return a list with the x position of each SiPM (in order).'''
        return [ sipm.X() for sipm in self.sipms ]

    def GetAllY( self ):
        ''' Return a list with the y position of each SiPM (in order).'''
        return [ sipm.Y() for sipm in self.sipms ]

    def OrderSiPMs(self):
        ''' Rearrange SiPMs by signal.'''
        self.sipms = OrderLike( self.GetAllQ(), self.sipms, 1 )
        return None

    def RemoveDiffuseLight( self, Type = 'absolute' ):
        ''' Erase SiPMs with less than QSIPMCUT.'''
        if Type == 'relative':
            qsipmcut = QSIPMCUTREL * self.Qmax()
        elif Type == 'absolute':
            qsipmcut = QSIPMCUTABS
        elif Type == 'both':
            self.RemoveDiffuseLight( 'absolute' )
            self.RemoveDiffuseLight( 'relative' )
            return None
        else:
            qsipmcut = 0.
        
        for i in Reverse( range( self.nsipms ) ):
            if Type=='absolute':
                self.sipms[i].Q -= qsipmcut
                if self.sipms[i].Q< 0.:
                    del self.sipms[i]
            elif Type=='relative':
                if self.sipms[i].Q < qsipmcut:
                    del self.sipms[i]
        self.nsipms = len( self.sipms )

    def ComputeCharge( self ):
        ''' Recalculate charge.'''
        self.Qcor = sum( self.GetAllQ() )

    def ComputeBaricenter( self, cut = 0. ):
        ''' Recalculate baricenter.'''
        if not self.nsipms:
            self.XB = self.YB = -100.
        
        X, Y, Q = self.GetAllX(), self.GetAllY(), self.GetAllQ()
        X, Y, Q = zip( *filter( lambda x: x[2] >= cut * Q[0], zip(X,Y,Q) ) )
        
        self.XB      = sum( [ x*q for x,q in zip(X,Q) ] )/sum(Q)
        self.YB      = sum( [ y*q for y,q in zip(Y,Q) ] )/sum(Q)
        self.RB      = ( self.XB**2 + self.YB**2 )**.5
        self.phiB    = Angle( self.XB, self.YB )
        self.coronaB = Corona( self.RB, R0corona )
        self.sectorB = PhiSector( self.phiB, Phi0corona )

    def ComputeRMS( self ):
        ''' Compute the RMS.'''
        self.RMSX = RMS( self.GetAllX(), self.GetAllQ(), self.X ) if self.nsipms > 1 else 0.
        self.RMSY = RMS( self.GetAllY(), self.GetAllQ(), self.Y ) if self.nsipms > 1 else 0.
        self.RMS  = sqrt( self.RMSX**2 + self.RMSY**2 )
    
    def ComputeCovarianceMatrix( self ):
        ''' Return the covariance matrix diagonalized: 1s matrix contains the eigenvalues and the second one the eigenvectors (in columns) .'''
        C = zeros(2,2)

        if self.nsipms < 4:
            print 'Not enough SiPMs to compute covariance matrix'
            self.CovarianceMatrix = None
            self.evals = None
            self.evecs = None
            self.RMS1  = None
            self.RMS2  = None
            self.VEC1  = None
            self.VEC2  = None
            self.theta = None
        
        X, Y, Q = self.GetAllX(), self.GetAllY(), self.GetAllQ()
        
        C[0][0] = Variance( X, Q, self.X )
        C[0][1] = Covariance( X, Y, Q, self.X, self.Y )
        C[1][0] = C[0][1]
        C[1][1] = Variance( Y, Q, self.Y )
        
        self.CovarianceMatrix = C
        evals, evecs = Diagonalize( C )
        self.RMS1    = math.sqrt( evals[0][0] )
        self.RMS2    = math.sqrt( evals[1][1] )
        self.VEC1    = evecs[0][0], evecs[1][0]
        self.VEC2    = evecs[0][1], evecs[1][1]
        self.theta   = atan( evecs[1][0] / evecs[0][0] ) * 180. / pi

    def ComputeSiPMsMatrix( self ):
        ''' Return the data of the SiPMs in matrix form.'''
        M = zeros( 151, 151 )
        for sipm in self.sipms:
            M[ int(sipm.X()) + 75 ][ int(sipm.Y()) + 75 ] = sipm.Q
        
        self.SiPMsMatrix = M
    
    def ComputeSource( self ):
        ''' Return the reconstructed source using the Fourier method.'''
        self.ComputeSiPMsMatrix()
        I = self.SiPMsMatrix
        S = IFT( Multiply2D( FT(I), GPSFT ) )
        self.source = Scale2D( S, self.E / Sum2D(S, absval = True) )
    
    def ComputeEcorrected( self, x = None, y = None ):
        ''' Return the corrected energy in coronas.'''
        self.Ecor = sum( [ pmt.Eweighted() for pmt in self.pmts.values() ] ) * CoronaCorrection[ self.corona ]

    def __repr__( self ):
        return '  E: {0}   Q: {1}   nsipms: {2}\n  SiPMs:\n{3}\n  PMTs:\n{4}'.format( self.E, self.Q, self.nsipms, '\n'.join( map( str, self.sipms ) ), '\n'.join( map( str, self.pmts.values() ) ) )

#SiPM instance
class sipm:
    def __init__( self, ID, Q, corona ):
        self.Q      = Q
        self.ID     = ID
        self.corona = corona;
    
    def X( self ):
        ''' X of the SiPM.'''
        return SiPM_map[self.ID][0]

    def Y( self ):
        ''' Y of the SiPM.'''
        return SiPM_map[self.ID][1]

    def R( self ):
        ''' R of the SiPM.'''
        return sqrt( self.X()**2 + self.Y()**2 )
    
    def distance( self, *XY ):
        ''' Return the distance of a SIPM relative to a point or to another SiPM.'''
        if isinstance( XY[0], (sipm, pmt) ):
            return self.distance( XY[0].X(), XY[0].Y() )
        
        return sqrt( ( self.X() - XY[0] )**2 + ( self.Y() - XY[1] )**2 )

    def __repr__( self ):
        return '    ID: {0}   Q: {1}   X: {2}  Y: {3}'.format( self.ID, self.Q, self.X(), self.Y() )


#PMT instance
class pmt:
    def __init__( self, ID, E ):
        self.ID = ID
        self.E  = E
    
    def Eweighted( self ):
        return self.E * PMT_weights[ self.ID ]
    
    def X( self ):
        ''' X of the PMT.'''
        return PMT_map[self.ID][0]

    def Y( self ):
        ''' Y of the PMT.'''
        return PMT_map[self.ID][1]

    def R( self ):
        ''' R of the PMT.'''
        return sqrt( self.X()**2 + self.Y()**2 )
    
    def distance( self, *XY ):
        ''' Return the distance of a SIPM relative to a point or to another SiPM.'''
        if isinstance( XY[0], (sipm, pmt) ):
            return self.distance( XY[0].X(), XY[0].Y() )
        
        return sqrt( ( self.X() - XY[0] )**2 + ( self.Y() - XY[1] )**2 )


    def __repr__(self):
        return '    ID: {0}   E: {1}   X: {2}  Y: {3}'.format( self.ID, self.E, self.X(), self.Y() )


#Data instance
class TrackBuilder(IReader):
    ''' This class is initialized with a list of names of files. It takes every file and opens it whenever necessary.'''
    def __init__(self, fname):
        IReader.__init__(self,fname)
        self.fname      = fname
        self.file       = None
        self.tracks     = None
        self.slices     = None
        self.sipms      = None
        self.pmts       = None
        self.Tcounter   = None
        self.Scounter   = None
        self.Sicounter  = None
        self.PMTcounter = None
        self.nevents    = None
    
    def open(self):
        self.file = ROOT.TFile(self.fname)

        # Reset variables.
        self.Reset()

    def close(self):
        self.file.Close()
    
    def eof(self):
        return self.Tcounter==self.nevents

    #Get a new track
    def read(self):
        
        # Move to the next file if the end is reached
        if self.Tcounter==self.nevents:
            self.eof()
    
    
        # Get an entry
        self.tracks.GetEntry(self.Tcounter)
        
        # Track variables, but sipms and pmts
        t = track(self.tracks.E,
                  self.tracks.Q,
                  self.tracks.nslices,
                  self.tracks.Emaxpos,
                  self.tracks.Qmaxpos,
                  self.tracks.isright,
                  [])

        for i in range(t.nslices):
            self.slices.GetEntry(self.Scounter+i)
            
            # Fill slice with its data
            s = slice(self.slices.E,
                      self.slices.Q,
                      self.slices.X,
                      self.slices.Y,
                      self.slices.Z,
                      self.slices.IDmax,
                      self.slices.nsipms,
                      [],{})
            
            # Fill sipms
            for j in range(s.nsipms):
                self.sipms.GetEntry(self.Sicounter+j)
                
                s.sipms += [ sipm(self.sipms.ID,
                                  self.sipms.Q,
                                  self.sipms.corona)]
            
            # Fill PMTs
            for j in range(NPMTS):
                self.pmts.GetEntry(self.PMTcounter+j)
                
                s.pmts[self.pmts.ID] = pmt( self.pmts.ID,
                                            self.pmts.E )
            
            # Order SiPMs and PMTs by charge-energy (decreasing)
            if s.nsipms:
                s.OrderSiPMs()
#                s.ComputeBaricenter()
#                s.ComputeRMS()
#                s.ComputeCovarianceMatrix()
#                s.ComputeSource()
#                s.ComputeEcorrected()

                
            # Update counters and add slice to the track.
            self.Sicounter  += s.nsipms
            self.PMTcounter += NPMTS
            t.slices        += [s]

        self.Scounter += t.nslices
        self.Tcounter += 1

        return t
    
    # Variables updater.
    def Reset(self):
        self.tracks     = self.file.Get('treetracks')
        self.slices     = self.file.Get('treeslices')
        self.sipms      = self.file.Get('treesipms' )
        self.pmts       = self.file.Get('treepmts'  )
        self.Tcounter   = 0
        self.Scounter   = 0
        self.Sicounter  = 0
        self.PMTcounter = 0
        self.nevents    = self.tracks.GetEntries()
        return None

    #Calculate the number of events 
    def totalevents(self):
        files = map( lambda x: ROOT.TFile(x), self.filelist )
        total = sum( map( lambda x: x.Get('treetracks').GetEntries(), files ) )
        map( lambda x: x.Close(), files )
        return total








