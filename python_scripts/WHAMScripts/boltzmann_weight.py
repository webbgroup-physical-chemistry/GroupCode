#!/usr/bin/python

nbins = 72

from numpy import array, sqrt, zeros, loadtxt
import sys
    
try :
    probfile = sys.argv[1]
    probabilities = loadtxt( probfile, "float64" )
    datafile = sys.argv[2]
    if sum(probabilities) > 1 or sum(probabilities) < 0.999 :
        print "\nWARNING\n>>>>Probabilities sum to %.30f<<<<\n" %sum(probabilities)
except:
    print "usage: ./script <prob file> <data file>"
    print "The probability file is the one-dimensional list of probabilities."
    print "The data file lists files:"
    print " <bin file> <data file>"
    print "which will be used for additional weighting of the result."
    sys.exit()

# the value of the degree of freedom in bin i which was sampled ncounts times is
# prob(i) * ( d[1] + d[2] + ... ) / ncounts

file = open( datafile )
datafilelines = file.readlines()
file.close()

binfile0, dofdata0 = datafilelines[0].split()
try :
    doffile0 = open( dofdata0 )
    dof0=doffile0.readlines()[0]
    values = len(dof0.split())
except IOError :
    print "file '%s' or '%s' not read, skipping" % ( binfile, dofdata )
    values = 1
grammar = "s" if values > 1 else ""
print "Looking at %i column%s of data" %(values,grammar)

for value in range(values) :
counts = [ [] for i in range(nbins) ]
nread = 0
    for line in datafilelines :
        binfile, dofdata = line.split()
        print "reading", dofdata, "using bins", binfile, ";",

        try : 
            binsfile = open( binfile )
            bins = binsfile.readlines()
            doffile = open( dofdata )
            dof=doffile.readlines()
    
            nread += 1
        except IOError :
            print "file '%s' or '%s' not read, skipping" % ( binfile, dofdata )
            continue

        lenbins = len(bins)
        lendof = len(dof)
        if lenbins != lendof :
            print "WARNING! %d points in bin file, %d points in dof file" % ( lenbins, lendof )
        else :
            print "%d points" % lenbins


        for i in range(len(bins)) :
            try :
                thisbin = int( bins[i].strip() )
                thisdof = float( dof[i].split()[value] )
            except ValueError :
                print "Incomprehensible '", dof[i].split()[value], "' at point", i
                continue
            except IndexError:
                print "IndexError in dof file %s at line %d" % ( dofdata, i )
                continue

            try:
                counts[thisbin].append( thisdof )
            except IndexError :
                print "IndexError"
                counts[0].append( thisdof ) # close enough

    total = 0
    vartotal = 0
    for i in range( nbins ) :
        #print "bin", i, "counts", len(counts[i])
        if len(counts[i]) > 0 :
            mean = array(counts[i]).mean()
            meansq = ( array( counts[i] )**2 ).mean()
            total += probabilities[i] * mean 
            vartotal += probabilities[i] * meansq

    print value,":",probfile, nread, total, sqrt(vartotal - total**2 )
