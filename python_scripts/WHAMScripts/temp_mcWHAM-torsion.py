#!/usr/bin/python

from Forces import *
from WHAMfunctions import *
from Distributions import WHAMMaxProposalDensity, montecarlo
from mpmath import log
from BefeusTools import TorsionExperiment
from numpy import array, dot, ones, savetxt, arange, inf
#from numpy import log
from pylab import figure, plot, show, bar, errorbar, axis, savefig
from cPickle import Pickler
import sys

#### MAIN
if __name__=="__main__":
    try : 
        filelistfile = open(sys.argv[1])
        begin, end, step, temperature = sys.argv[2:6]
    except :
        print "usage: ./script <file name> <begin> <end> <step> <temperature>"
        sys.exit()

    try :
        outputdir = sys.argv[6]
    except:
        print "no output directory supplied, using '.'"
        outputdir = "."

    name = sys.argv[1]

    begin = float( begin )
    end = float( end )
    step = float( step )
    kT = float( .00831447 )*float( temperature )
    beta = 1/float( kT )

    experiments=[ TorsionExperiment(line) for line in filelistfile.readlines() if line[0] != "#" ]
    filelistfile.close()

    nexperiments = len(experiments)

    for experiment in experiments : 
        experiment.bin( begin, end, step )
        print "*** finished reading", experiment
        print " counts =", experiment.counts

    omegas = array( [ experiment.omegas for experiment in experiments ] )
    print "Built (%dx%d) omega matrix" % omegas.shape

    # count the state counts
    # this may have to change if there are states without counts (esp
    # given a Jeffreys prior)
    nstates = len( omegas[0] )
    counts = zeros( nstates, "float" )
    samples = zeros( nexperiments, "float" )
    for i in range( nexperiments ):
        experiment = experiments[i]
        counts += experiment.counts
        samples[i] = experiment.counts.sum()
    print "State counts:", counts
    print "Samples:", samples

    # settings and initial conditions
    w0 = ones( nstates, "float" )/nstates
    force = WHAMForce( counts, omegas, samples )
    density = WHAMDensity( counts, omegas, samples )

    # optimize
    niter = 50 
    tol = 10**-9 
    print "OPTIMIZING"
    opttrajectory = DoWHAM(w0, counts, samples, omegas, niter, tol )
    print "Ran %d points of optimization" % len(opttrajectory)
    print opttrajectory[-1]

    # sample
    niter = 10 
    diriscale = 10**8
    print "SAMPLING"
    targetDensity = WHAMDensity( omegas = omegas, counts = counts, samples = samples )
    proposalDensity = WHAMMaxProposalDensity( opttrajectory[-1].tolist(), diriscale )
    # make empty bins useful. can't do this with proposalDensity
    counts += ones(nstates, "float")
    trajectory = montecarlo( opttrajectory[-1], proposalDensity, targetDensity, niter ) 

    # results
    #print "last point =", trajectory[-1]
    print len(trajectory), "points sampled"
    avg = trajectory[1:].mean(axis=0)
    std = trajectory[1:].std(axis=0)
    #print "average =",  avg
    #print "standard deviation =", std
    
    # calculate the mean of -kT ln(w) and its second moment
    # Rule: for any sample that had '0', assume that the mean
    # is 0 and the variance is infinite.
    mean_potential = zeros( nstates, "float" )
    std_potential = zeros( nstates, "float" )
    trajectoryT = trajectory.transpose()
    for state in range( nstates ):
        statesamples = trajectoryT[ state ]
        if any( statesamples == 0 ):
            mean = 1000 # None
            std = inf#  None 
        else :
            tot = mpf(0.0)
            tot2 = mpf(0.0)
            for sample in statesamples :
                logsample = -kT*log(sample)
                tot += logsample
                tot2 += logsample*logsample
            
            mean = tot/niter 
            mom2 = tot2/niter
            std =  sqrt( mom2 - mean*mean )

        mean_potential[state] = mean
        std_potential[state] = std.real # danger!

    xvalues = arange( begin, end, step ) 

    # display
    figure(1)
    errorbar( xvalues, mean_potential, yerr=2*std_potential )   
    axis([ begin,end,0, mean_potential.max()*1.1] )
    savefig( "%s/%s.eps" % ( outputdir, name) ) 

    figure(2)
    plot( trajectory )

    #show()

    # save
    savetxt( "%s/%s.traj" % ( outputdir, name ), trajectory )
    savetxt( "%s/%s.mean" % ( outputdir, name ), mean_potential )
    savetxt( "%s/%s.std" % ( outputdir, name ) , std_potential )
