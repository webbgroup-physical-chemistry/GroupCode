import unittest
from mpmath import exp, mpf
from numpy import array, dot, ones, inf, zeros
from numpy.random import normal, shuffle, random, exponential, randint
from numpy.random.mtrand import dirichlet, multinomial
import sys
from time import time
## TARGET DISTRIBUTIONS
# Parent class
class TargetDistribution :
    def __add__( self, other ): return Mixture( self, other )

# GAUSSIAN
class Gaussian( TargetDistribution ) :
    """A class for a Gaussian target density."""
    def __init__( self, mu, sigma ):
        self.mu = mu
        self.sigma = sigma
        self.invsigma = 1.0/sigma
        self.sigma2 = sigma*sigma
        self.halfinvsigma2 = 0.5/self.sigma2

    def __call__(self, x):
        return self.invsigma*exp( -(x-self.mu)*(x-self.mu)*self.halfinvsigma2 )

class Mixture( TargetDistribution ) :
    """A class to deal with sums of distributions.""" 
    def __init__( self, distribution1, distribution2 ):
        self.distributions = ( distribution1, distribution2 )
    
    def __call__( self, x ):
        results = [ dist(x) for dist in self.distributions ]
        results = array(results)
        return results.sum()

# EXPONENTIAL
class Exponential( TargetDistribution ) :
    """Exponential target density."""
    def __init__( self, tau ):
        self.tau = tau  
        self.neginvtau = -1.0/tau
        
    def __call__( self, x ):
        return exp( self.neginvtau*x )

# WHAM Density
class WHAMDensity( TargetDistribution ):
    """The posterior distribution from a WHAM calculation. Three numpy arrays are needed:\n\tomegas - the matrix of integrated Boltzmann weights. These should probably be mpf types!\n\tcounts - the total number of counts of each state in all simulations.\n\tsamples - the total number of samples in each simulation."""
    def __init__( self, omegas, counts, samples ):
        self.omegas = array(omegas, mpf)
        self.counts = array(counts)
        self.samples = array(samples)

        self.nexperiments, self.nstates = self.omegas.shape

        #self.omegaColumns = array( [ self.omegas[:,i] for i in range(self.nstates) ] )
        self.firstcolumn = self.omegas[:,0]

        # this is the 1st row minus the mth row, from 2 .. nexperiments
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

        # Numpy version
        numerator = (wNotFirst**self.countsNotFirst).prod()
        numerator *= ( 1 - wNotFirst.sum() )**self.countsFirst

        # Alternate (loop) version
        """
        numerator = ( 1 - wNotFirst.sum() )**self.counts[0]
        for a in range(1,self.nstates ) :
            numerator *= w[a]**self.counts[a]
        """
        # Denominator: Numpy version
        # this is slow though too because of mpf!
        """
        denominatorFactors = array([ mpf(dot(row[1:],wNotFirst)) for row in self.deltaOmegas ])
        denominatorFactors += self.firstcolumn
        denominatorFactors = denominatorFactors**self.samples
        denominator = denominatorFactors.prod()
        """
        # Alternate (loop) version
        denominator = 1
        for m in range(self.nexperiments ):
            total = self.omegas[m,0]
            for a in range(1,self.nstates ):
                deltaOmega = self.omegas[m,a] - self.omegas[m,0]
                total += mpf(deltaOmega)*w[a]
            denominator *= total**self.samples[m]

        return numerator/denominator
    
    # really clunky
    def optimize( self, x0, nsteps ):
        dim = len(x0)
        proposalDensity = RandomDirichletProposalDensity( dim )
        #proposalDensity = DirichletSwapProposalDensity( dim )
        trajectory = zeros( (nsteps,dim), "float" )
        trajectory[0] = x0
        print "Optimizing (start at",x0,") ..."
        i = 0
        count = 0
        while True :
            try :
                oldpoint = trajectory[i]
                #proposal = proposalDensity( oldpoint )
                proposal = proposalDensity()
                if self(proposal) > self(oldpoint):
                    trajectory[i+1] = proposal
                    i += 1
                count += 1
                if count % 1000 == 0 : print "step", count, i, "accepted, current state", trajectory[i-1]
            except IndexError :
                break
        return trajectory

class WHAMUnitTest( unittest.TestCase ):
    def setUp( self ):
        self.nexperiments = 5
        self.nstates = 10

    def setup_ones( self, scale1 = 1, scale2 = 1, scale3 = 1   ):
        omegas = scale1*ones( (self.nexperiments,self.nstates), "float" )
        counts = scale2*ones( self.nstates, "float" )
        samples= scale3*ones( self.nexperiments, "float" )
        self.distribution = WHAMDensity( omegas, counts, samples )

    """
    def testcall_zeros( self ):
        self.setup_ones( scale1 = 0 )
        data = ones( self.nstates, "float")/self.nstates
        # for some reason this actually returns inf
        #self.assertRaises(ZeroDivisionError, self.distribution, data )
        self.assertEqual( inf, self.distribution(data) )
    """
    
    def testcall_ones( self ):
        self.setup_ones( )
        data = ones( self.nstates, "float")/self.nstates
        expected = self.nstates**-self.nstates
        self.assertAlmostEqual( expected, self.distribution(data) )

    def testcall_onetenth( self ):
        self.setup_ones( scale1 = 0.1 )
        data = ones( self.nstates, "float" )/self.nstates 
        expected = self.nstates**-self.nstates * 0.1**-self.nexperiments
        #print expected
        #print self.distribution(data)
        self.assertAlmostEqual( expected, self.distribution(data) )

# PROPOSAL DENSITIES 
# A symmetric Gaussian proposal density with user-defined width.
# This is probably slow.
class GaussianProposalDensity :
    def __init__( self, param ): self.param = param
    def __call__( self, x ): return normal(x,self.param)

# A Gaussian proposal density for x >= 0.
class GaussianPositiveProposalDensity :
    def __init__( self, param ): self.param = param
    def __call__( self, x ):
        proposal = -1
        while proposal < 0 : proposal = normal(x,self.param)
        return proposal

# Dirichlet proposal density to give a vector
# of length M whose elements sum to 1. The vector
# is chosen completely at random.
class RandomDirichletProposalDensity :
    def __init__( self, M ) :
        self.M = M
        self.alphas = ones( M )
    def __call__( self, nothing=None ): return dirichlet( self.alphas )

class RandomDirichletProposalUnitTest( unittest.TestCase ):
    def setUp( self ):
        M = 10
        self.proposalDensity = RandomDirichletProposalDensity(M)

    def test_rdp( self ):
        for i in range( 10 ):
            self.assertAlmostEqual( 1, self.proposalDensity().sum() )

# Pick two elements at random and let one be x and the other be 1-x
class DirichletSwapProposalDensity :
    def __init__( self, M ):
        self.indices = range(M)

    def __call__( self, w ) :
        shuffle( self.indices )
        index1, index2 = self.indices[:2]
        newr1 = random()
        newr2 = 1-newr1
        neww = w.copy()
        neww[index1] = newr1
        neww[index2] = newr2
        return neww
        
class DirichletSwapProposalDensityUnitTest( unittest.TestCase ) :
    def setUp( self ):
        self.M = 10
        self.proposalDensity = DirichletSwapProposalDensity(self.M)

    # this one should succeed if w and proposal are different
    def test_get_w_back( self ):
        w = ones( self.M, "float" )
        proposal = self.proposalDensity(w)
        #print w
        #print proposal
        #print (w==proposal).all()
        self.assertFalse( (w==proposal).all() )

    def examples( self ):
        w = ones( self.M, "float")/self.M
        for i in range(10 ):
            w = self.proposalDensity(w)
            print w

# Change half of states by (r), half by -r
# If odd, just leave one unchanged.
class ShiftProposalDensity :
    def __init__( self, k, nstates ):
        self.k = k
        self.nstates = nstates

        # make an array with half the entries 1, half -1, and any remaining zero
        half = nstates/2
        myones = ones( half, "int" )
        mynegs = -ones( half, "int" )
        values = zeros( nstates, "int" )
        values[0:half] = myones
        values[half:2*half] = mynegs
        self.values = values

    def __call__( self, w ) :
        shuffle(self.values)
        move= exponential(self.k) * self.values
        neww = w + move
        while not self.validrange( neww ) :
            shuffle(self.values)
            move= exponential(self.k) * self.values
            neww = w + move
        return neww

    def validrange( self, arr ):
        validity = ( arr >=0 ) * ( arr <= 1 )
        return validity.all()

class ShiftProposalDensityUnitTest( unittest.TestCase ) :
    def setUp( self ):
        self.k = .1 
        self.nstates = 10
        self.proposalDensity = ShiftProposalDensity(k=self.k, nstates=self.nstates)
        self.warray1 = ones(self.nstates, "float") /self.nstates

        self.warray2 = self.warray1.copy()
        self.warray2[0] = self.warray1[0]*.5
        self.warray2[1] = self.warray1[1] + self.warray1[0]*.5

    def validrange( self, arr ):
        validity = ( arr >=0 ) * ( arr <= 1 )
        return validity.all()

    # make sure the test arrays make sense
    # both == 1
    # 0 <= both <= 1
    def test_sums( self ):
        self.assertAlmostEqual( 1, self.warray1.sum() )
        self.assertAlmostEqual( 1, self.warray2.sum() )

    def test_ranges( self ):
        self.assertTrue( self.validrange( self.warray1 ) )
        self.assertTrue( self.validrange( self.warray2 ) )
        
    # now test whether the proposal density yields a valid distribution
    def test_proposal_sum( self ):
        proposal1 = self.proposalDensity( self.warray1 )
        self.assertAlmostEqual( 1, proposal1.sum() )
        proposal2 = self.proposalDensity( self.warray2 )
        self.assertAlmostEqual( 1, proposal2.sum() )

    # test whetehr the proposal density yields values 0<=x<=1
    def test_proposal_valid( self ):
        for i in range(100) :
            proposal1 = self.proposalDensity( self.warray1 )
            self.assertTrue( self.validrange( proposal1 ) )
            proposal2 = self.proposalDensity( self.warray2 )
            self.assertTrue( self.validrange( proposal2 ) )

# Use the maximum likelihood estimate (regular WHAM) probability vector
# as defining a multinomial distribution. Sample K times from that 
# multinomial. The proposal density is n_i/K. The advantage is that this
# proposal density is independent of the current step, giving an 
# independence chain. 
class WHAMMaxProposalDensity :
    def __init__( self, pvector, Msamples ):
        self.pvector = array(pvector)
        self.Msamples = Msamples
        self.Msamplesfloat = float(self.Msamples)

    def type( self ): return self.pvector.dtype

    def __call__( self, nothing=None ) :
        nsample = randint( 1, self.Msamples )
        fnsample = float(nsample)
        sample = multinomial( nsample, self.pvector )/fnsample
        #print sample
        return sample

class test_WHAMMaxProposalDensity( unittest.TestCase ):
    def setUp( self ):
        self.density = WHAMMaxProposalDensity( (.1,.3,.6), 10 )

    def test_type( self ): print self.density.type()
    
    def test_sample( self  ):
        for i in range(10):
            print self.density()

def montecarlo( w0, proposalDensity, targetDistribution, niter ) :
    # Initialize the trajectory
    nstates = len(w0)
    trajectory = zeros( (niter,nstates), "float" )
    trajectory[0] = w0

    # sample
    print "sampling %d points ..." % niter
    i = 0
    naccept = 0
    t0 = time()
    while True :
        try: 
            # old x 
            oldx = trajectory[i] 

            # propose a new x
            newx = proposalDensity( oldx )

            oldweight = targetDistribution( oldx )
            newweight = targetDistribution( newx )

            # compute the transition probability
            rcut = newweight/oldweight

            # update
            #print i, newweight,"/",oldweight,"=",rcut
            if random() < rcut: 
                #print "step", i, "accepted with prob", rcut
                trajectory[i+1] = newx
                naccept += 1
            else:
                #print "step", i, "rejected with prob", rcut
                trajectory[i+1] = oldx
            i += 1
        except IndexError :
            print "sampling ended after the %dth point" % i
            print "accepted %d of %d points (%.3f)" % ( naccept, niter, float(naccept)/niter )
            break

    return trajectory

# For unittesting
if __name__ == "__main__": unittest.main()
