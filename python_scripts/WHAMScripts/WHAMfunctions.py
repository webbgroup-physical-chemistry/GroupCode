import unittest
from mpmath import mpf, mpmathify
from numpy import array, dot, ones, zeros
import sys

def buildTest():
    nstates = 4
    nsimulations = 10
    omegas = ones( (nsimulations, nstates), "float" )
    wvec = ones( nstates, "float")/nstates
    Nm = 10*ones( nsimulations, "float" )
    return omegas, Nm, wvec

def fm( wvec, omegarow, Nm ): return Nm/dot( omegarow, wvec )

def WHAMstep( wvec, ncounts, nobs, omegas ):
    nstates = len(wvec)
    nexperiments = len(nobs)

    fms = zeros( nexperiments, "float" )
    for m in range(nexperiments) : fms[m] = fm( wvec, omegas[m], nobs[m] ) 

    wguess = zeros( nstates, "float" )
    for a in range(nstates):
        omegacol = omegas[:,a]
        fattractive = mpf( dot( omegacol, fms ) )
        wguess[a] = mpf(ncounts[a])/fattractive
  #      print a,wguess[a],"=",ncounts[a],"/",fattractive,"=",omegacol,"dot",fms

    wguess/=wguess.sum() 

    return wguess

def DoWHAM( wvec, ncounts, nobs, omegas, niter, tol=10**-6 ):

    # do one wham step
    trajectory = zeros( (niter,len(wvec)), mpf )
    trajectory[0] = wvec
    i = 0
    while True :
        try :
            trajectory[i+1] = WHAMstep( trajectory[i], ncounts, nobs, omegas )

            # compute % change and quit if acceptable
            # should this be mean() instead of sum()
            if (( trajectory[i+1]-trajectory[i] )**2).sum() < tol :
                print "converged to", tol,"at", i
                trajectory = trajectory[:i+1]
                break

            if 0/trajectory[i+1].sum() != 0 :
                print "nan at", i, trajectory[i]
                trajectory = trajectory[:i+1]
                break
            i += 1
        except IndexError :
            break

    return trajectory

# convert a numpy array to mp type
def make_mp( arr ) :
    return array( [ [ mpmathify( row[i] ) for i in range(len(row)) ] for row in arr ] )

class Test_make_mp( unittest.TestCase ):
    def setUp( self ):
        self.arr = ones( ( 3,5), "int" )

    def test_make_mp( self ):
        mparr = make_mp( self.arr )
        self.assertEquals(  mparr.dtype, object )

if __name__=="__main__":
    unittest.main()
