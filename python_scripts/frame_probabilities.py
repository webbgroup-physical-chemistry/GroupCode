#!/usr/bin/python

#nbins = 72

from numpy import array, sqrt, zeros, loadtxt
import sys
    
try :
    probfile = sys.argv[1]
    probabilities = loadtxt( probfile, "float64" )
    nbins = len(probabilities)
    datafile = sys.argv[2]
    if sum(probabilities) > 1 or sum(probabilities) < 0.9999 :
        print "\nWARNING\n>>>>Probabilities sum to %.10f<<<<\n" %sum(probabilities)
except:
    print "usage: ./script <prob file> <data file>"
    print "The probability file is the one-dimensional list of probabilities."
    print "The data file lists files:"
    print " <bin file> <data file>"
    print "which will not use the <data file>, it is just kept to keep"
    print "input consistant with boltzmann_weight.py."
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
    print "file '%s' or '%s' not read, skipping" % ( binfile0, dofdata0 )
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
    
            nread += 1
        except IOError :
            print "file '%s' or '%s' not read, skipping" % ( binfile, dofdata )
            continue

        lenbins = len(bins)
        print "%d points" % lenbins


        for i in range(len(bins)) :
            try :
                thisbin = int( bins[i].strip() )
            except :
                print "cannot read %s, exiting..." %bins
                sys.exit()

            try:
                counts[thisbin].append( '0' )
            except IndexError :
                print "IndexError"
                counts[0].append( '0' ) # close enough

    

binprob=[ [] for i in range(nbins) ]
for i in range( nbins ) :
    #print "bin", i, "counts", len(counts[i])
    if len(counts[i]) > 0 :
        binprob[i] = probabilities[i] / len(counts[i])
    else :
        binprob[i] = probabilities[i]



n2chi={}
z=0
for r1 in xrange(0,360,30) :
    for r2 in xrange(0,360,30) :
        n2chi[z]="%s-%s" %(r1,r2)
        z+=1

g=[]
for line in datafilelines :
    binfile,dofdata = line.split()
    binname = binfile.split('.')
    if binname[1] == "" : x="%s..%s"%(binname[0],binname[2])
    else : x=binname[0]
    probfilename = "%s_%s.prob" %( x, n2chi[int(binname[-2])] )
    probfile = open(probfilename,"w")
        
    try :
        binsfile=open(binfile)
        bins=binsfile.readlines()
        binsfile.close()
    except :
        print "file '%s' not read, exiting" %binfile
        sys.exit()
    for line in bins :
        this_bin=int(line.strip())
        probfile.write("%.9f\n" %binprob[this_bin])
        g.append("%.9f" %binprob[this_bin])
    probfile.close()


