#!/usr/bin/python

from mpmath import mpf, sqrt, power
from numpy import abs, array, dot, inf, ones, zeros
from numpy.random import normal
import unittest
from time import time

mpf1= mpf(1.0)

## The force on a WHAM particle.
class WHAMForce :
    def __init__( self, counts, omegas, samples ):
        self.counts = counts
        self.nstates = len(counts)
        self.omegas = omegas
        self.nsimulations = len(omegas)
        self.nsamples = samples

        # compute deltaOmegas
        # compute omegaColumns
        self.omegaColumns = [ omegas[:,i] for i in range( omegas.shape[1] ) ] # correspond to states
        self.deltaOmegas = self.omegaColumns-self.omegaColumns[0]

    def __call__( self, w ):
        counts = self.counts
        nstates = self.nstates
        deltaOmegas = self.deltaOmegas
        omegaColumns = self.omegaColumns
        nsimulations = self.nsimulations
        nsamples = self.nsamples

        # functions used by the attractive and interactive parts
        wsubset = w[1:]
        ffunction = 1- wsubset.sum()

        # Repulsive parts n_a/w_a
        repulsivePart = counts/w 

        # Attraction to normalization
        attractivePart = -float(counts[0] / ffunction )
    
        # Interaction between weights
        interactionPart = zeros(nstates, "float")
        dotprod = [ dot( deltaOmegas[:,i][1:], wsubset ) for i in range(nsimulations)]
        denominator = array(omegaColumns[0] + dotprod, "float" )
        factor = nsamples/denominator 
        interactionPart = -array([ dot( row, factor ) for row in deltaOmegas ], mpf )
        force = repulsivePart + attractivePart + interactionPart
        force[0] = 0
        return force

class WHAMForceUnitTest( unittest.TestCase ) :
    def setUp( self ):
        self.nsimulations = 10
        self.nstates = 5

    def setup_ones( self, scale1 = 1, scale2 = 1, scale3 = 1   ):
        omegas = scale1*ones( (self.nsimulations,self.nstates), "float" )
        counts = scale2*ones( self.nstates, "float" )
        samples= scale3*ones( self.nsimulations, "float" )
        self.force = WHAMForce( counts, omegas, samples )

    """
    def testcall_zeros( self ):
        self.setup_ones( scale1 = 0 )
        data = ones( self.nstates, "float")/self.nstates
        self.assertEqual( inf, self.force(data) )

    """
    def testcall_ones( self ):
        self.setup_ones( )
        data = ones( self.nstates, "float")/self.nstates
        #self.assertAlmostEqual( expected, self.force(data) )
        print data
        print self.force(data)

    """
    def testcall_onetenth( self ):
        self.setup_ones( scale1 = 0.1 )
        data = ones( self.nstates, "float" )/self.nstates
        expected = self.nstates**-self.nstates * 0.1**-self.nsimulations
        self.assertAlmostEqual( expected, self.force(data) )
    """

# A Brownian dynamics integrator.
# The 'temperature' might not be correct here (but it's close).
def integrate( w0, force, nsteps, tau, dim ) :
    sqrt2tau = sqrt( 2*tau )
    trajectory = zeros( (nsteps,dim), "float" )
    trajectory[0] = w0
    i = 1
    while 1 :
        try :
            currentPoint = trajectory[i-1]
            randomForce = normal(0,sqrt2tau,dim)
            systematicForce = tau * force( currentPoint )
            newPoint = currentPoint + randomForce + systematicForce
            newPoint[0] = 1 - newPoint[1:].sum()
            #trajectory[ i ] = newPoint 

            # A hack!
            newPoint = abs(newPoint)
            trajectory[ i] = newPoint/newPoint.sum()

            i += 1    
        except IndexError :
            if i == trajectory.shape[0] : break
            else : raise
    return trajectory

# An optimizer, just like the integrator above, but without a stochastic
# component. 
# To - do: test the density of the new point, and reject if it's less than
# the old point. Scale dt down and try again. 
def findmax( w0, force, density, nsteps, tau, tauscale, nstates, tol=10**-5, maxtauscale=100 ) :
    print "starting optimization at", w0
    ntauscale = 0
    trajectory = zeros( (nsteps,nstates), mpf )
    trajectory[0] = w0
    i = 1
    while 1 :
        accept = True
        try :
            currentPoint = trajectory[i-1]
            systematicForce = tau * tauscale**ntauscale * force( currentPoint )
            newPoint = currentPoint + systematicForce 
            newPoint[0] = 1 - newPoint[1:].sum()
           
            # check if this is a good point
            if (newPoint<0).any() or (newPoint>1).any() :
                print "Bad sample", newPoint[:10], "... selected, rescaling step size"
                accept = False
                ntauscale += 1           

            # reject and scale if this isn't a better point
            olddensity = density(currentPoint)
            newdensity = density(newPoint)
            if newdensity < olddensity :
                print "new density", newdensity, "< old density", olddensity, "-- not accepting, and rescaling step size"
                accept = False
                ntauscale += 1

            # quit if we've scaled tau too many times
            if accept and ntauscale >= maxtauscale :
                print "We have scaled tau", ntauscale, "times"
                print "Final tau scaling was", tauscale,"**", ntauscale, "after", i-1, "steps"
                print "Density change was", newdensity-olddensity, "from", olddensity, "to", newdensity
                break
        
            if accept :
                # check to see if we might have converged
                if abs(newdensity - olddensity)/olddensity < tol and newdensity > 0 :
                    print "Density change less than", tol, "after %d steps" % i
                    break
                trajectory[ i ] = newPoint 
                i += 1    

        except IndexError :
            if i == trajectory.shape[0] : break
            else : raise
    print "Optimization ended after", i,"steps"
    return trajectory[:i]

# A WHAM density for the optimizer findmax() above. The WHAM force is the derivative
# of the log of this density. It would be better to link the two more closely instead
# of two seperate classes. 
class WHAMDensity :
    """The posterior distribution from a WHAM calculation. Three numpy arrays are needed:\n\tomegas - the matrix of integrated Boltzmann weights. These should probably be mpf types!\n\tcounts - the total number of counts of each state in all simulations.\n\tsamples - the total number of samples in each simulation."""
    def __init__( self, counts, omegas, samples ):
        self.omegas = array(omegas, mpf)
        self.counts = array(counts)
        self.samples = array(samples)

        self.nexperiments, self.nstates = self.omegas.shape
        self.firstcolumn = self.omegas[:,0]
        firstrow = self.omegas[0]
        self.deltaOmegas = [ self.omegas[m] - firstrow for m in range(self.nexperiments) ]

        self.countsFirst = self.counts[0]
        self.countsNotFirst = self.counts[1:]
        self.samplesNotFirst = self.samples[1:]

    # This assumes a vector of length self.nsamples. Therefore, the first element
    # of this vector is not important, since it's fixed by the normalization condition
    # for w, ie, sum(w[m]) = 1. However, we do not check that the normalization is 
    # satisfied -- that should be the job of the proposal density. 
    def __call__(self, w ):
        wNotFirst = w[1:]

        numerator = array( [ power( wNotFirst[i], self.countsNotFirst[i] ) for i in range( self.nstates-1 ) ] ).prod()
        #numerator = (wNotFirst**self.countsNotFirst).prod()
        numerator *= ( mpf1 - wNotFirst.sum() )**self.countsFirst

        denominator = mpf1
        for m in range(self.nexperiments ):
            total = self.omegas[m,0]
            for a in range(1,self.nstates ):
                deltaOmega = self.omegas[m,a] - self.omegas[m,0]
                total += mpf(deltaOmega)*w[a]
            denominator *= total**self.samples[m]

        return numerator/denominator

if __name__ == "__main__" : unittest.main()
