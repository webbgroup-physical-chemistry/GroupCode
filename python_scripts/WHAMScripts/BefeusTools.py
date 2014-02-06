import unittest
from mpmath import pi, exp
from numpy import arange,array, histogram, zeros, arccos
from math import ceil, cos, sin
DEG2RAD = pi/180.0
    
# Added by Andrew 8/15/2012
def angle_difference( a, b):
    #"The arccos( dot(a,b) ) always returns the absolute value of the angle difference.  Useful here!"
    return arccos( cos(a*pi/180)*cos(b*pi/180)+sin(a*pi/180)*sin(b*pi/180))*180/pi

class Experiment :
    def readtrajectory( self ):
        try : file = open( self.filename )
        except IOError : 
            print "file %s not found" % self.filename
            return 
            
        data = []
        for line in file.readlines() :
            if line[0] != "#" :
                data.append( float( line.split()[-1] ) )
        file.close()
        self.trajectory = array(data)
        self.npoints = len(self.trajectory.data) # this is N_m

    def bin( self, begin, end, step ):
        self.bins = arange( begin, end+step, step )
        self.counts = histogram( self.trajectory, self.bins )[0]
        self.nbins = len(self.counts) # self.bins has one extra entry, because it defines 
                                      # the right-hand bin wall 
        # Determine which bins were visited. These are the only bins that this experiment
        # can illuminate, according to the Befeus procedure.
        self.binsVisited = []
        self.omegas = zeros( self.nbins, "object" ) 
        for i in range( self.nbins ) :
            if self.counts[i] > 0 : 
                self.binsVisited.append(i)
            a = i*step + begin
            b = (i+1)*step + begin
            self.omegas[i] = self.integrate(a,b)

class LinearExperiment( Experiment ) :
    def __init__( self, line ):
        filename, reference, forceconstant, beta = line.split()
        self.filename = filename
        self.reference = float(reference)
        self.forceconstant = float(forceconstant)
        self.beta = float(beta)
        self.halfforceconstant = 0.5*self.forceconstant
        self.readtrajectory()

    def __str__( self ):
        return self.filename

    # trapezoidal integration
    def integrate( self, a, b ):
        halfk = self.halfforceconstant
        x0 = self.reference
        beta = self.beta
        avalue = exp( -beta*halfk*(a-x0)**2 )
        bvalue = exp( -beta*halfk*(b-x0)**2 )
        return 0.5*(b-a)*( avalue + bvalue ) 

class TestLinearExperiment( unittest.TestCase ):
    def setUp( self ):
        self.experiment = LinearExperiment( "null 0 1 1" )
        
    def test_integrate( self ):
        self.assertAlmostEqual( self.experiment.integrate(0,1), 0.803265, places=6 )

### TORSION ###
class NewTorsionExperiment( Experiment ):
    def __init__( self, line ):
        filename, phi0, deltaphi, kfac = line.split()[:4]
        self.filename = filename
        self.phi0 = float(phi0)
        self.deltaphi = float(deltaphi)
        self.kfac = (pi/180)*(pi/180)*float(kfac)

class TorsionExperiment( Experiment ) :
    def __init__( self, line ) :
        filename, phi0, deltaphi, kfac, beta = line.split() # kludge city

        self.filename = filename
        self.phi0 = float(phi0)
        self.deltaphi = float(deltaphi)
        self.beta = float(beta)
        self.readtrajectory() 

        # convert force constant from energy/radian**2 to energy/degree**2
        self.kfac = float(kfac) * (pi/180) * (pi/180)
      
        print "working on", self.filename
 
    def __str__( self ):
        return "Experiment %s, phi0= %f, deltaphi= %f, kfac= %f, beta= %f" % ( self.filename, self.phi0, self.deltaphi, self.kfac, self.beta )

    # Do the integral from a to b of the GROMACS dihedral restraint function.
    # The Boltzmann-weighted potential is in self.weight_dihedral_restraint_potential().
    # This is a (naive?) trapezoidal integration. (Kincaid and Cheney p. 481)
    def integrate( self, a, b ) :
        function = self.weight_dihedral_restraint_potential
        #print "integrating from", a,"to", b,"degrees"
        #print a*DEG2RAD, "to", b*DEG2RAD, "radians"
        #print "function(a) = ",function(a),"and function(b) =", function(b)
        #return 0.5*(b-a)*DEG2RAD*( function(a) + function(b) )
        return 0.5*(b-a)*( function(a) + function(b) )

    def weight_dihedral_restraint_potential( self, phi ):
        phi0 = self.phi0
        deltaphi = self.deltaphi
        forceconst = self.kfac
##        diffphi = phi-phi0
 
        # This is copied from GROMACS src/gmxlib/dihres.c.
        # It's not quite right, since it only accounts for dihedrals as much as
        # 2 pi beyond the allowed range of (-pi,pi) (?).
        # Ignored because GMX outputs are always (-pi,pi).
##        if diffphi >= 180 : diffphi -= 360
##        elif diffphi < -180 : diffphi += 360
                
##        if diffphi > deltaphi :
##            ddp = diffphi - deltaphi 
##        elif diffphi < -deltaphi :
##            ddp = diffphi + deltaphi
##        else : ddp = 0
        
        diffphi = angle_difference(phi,phi0)
        if diffphi > deltaphi :
            ddp = angle_difference(diffphi, deltaphi)
        else :
            ddp=0

##        ddp *= DEG2RAD 
            
        energy = 0.5*forceconst*ddp*ddp
        return exp(-self.beta*energy)

    def dihedral_restraint_potential( self, phi ):
##        diffphi = phi-phi0
        
        # This is copied from GROMACS src/gmxlib/dihres.c.
        # It's not quite right, since it only accounts for dihedrals as much as
        # 2 pi beyond the allowed range of (-pi,pi) (?).
##        if diffphi >= 180 : diffphi -= 360
##        elif diffphi < -180 : diffphi += 360
        
##        if diffphi > deltaphi :
##            ddp = diffphi - deltaphi 
##        elif diffphi < -deltaphi :
##            ddp = diffphi + deltaphi
##        else : ddp = 0
        
        diffphi = angle_difference(phi,phi0)
        if diffphi > deltaphi :
            ddp = angle_difference(diffphi, deltaphi)
        else :
            ddp=0
        
##        ddp *= DEG2RAD 
        energy = 0.5*forceconst*ddp*ddp
        return energy
    


if __name__ == "__main__" :
    unittest.main() 


