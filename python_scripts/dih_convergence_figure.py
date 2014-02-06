#! /usr/bin/python

from pylab import *
from optparse import OptionParser
import sys

print "\n    This is used to make a .pdf figure comparing the 1st half and last half dihedral distributions to the total dihedral distribution.  It requires the existance of each <file.prob> prior to execution.\n\n"
parser=OptionParser()
parser.add_option("-a", "--full", dest="full", help="Full simulation probability file")
parser.add_option("-f", "--half1", dest="half1", help="First half of simulation probability file")
parser.add_option("-s", "--half2", dest="half2", help="Second half of simulation probability file")
parser.add_option("-o", "--out", dest="name", help="Output file pdf name", default="out")
parser.add_option("-x", "--chi", dest="chi", help="Chi <angle>", default=1)
parser.add_option("-t", "--title", dest="pdftitle", help="Figure title", default="")
parser.add_option("-g", "--grid", dest="grid", help="True: Enable axis labels.  Default: False", default=False)
(options,args)=parser.parse_args()

fullfile=options.full
half1=options.half1
half2=options.half2
name=options.name
chi=options.chi
pdftitle=options.pdftitle
grid=options.grid


try :
    full=open(fullfile)
    fulllines=full.readlines()
    full.close()
    h1=open(half1)
    h1lines=h1.readlines()
    h1.close()
    h2=open(half2)
    h2lines=h2.readlines()
    h2.close()
except :
    print "usage: %s -h" %sys.argv[0]
    sys.exit()

    
if len(h1lines) != len(h2lines) != len(fulllines) :
    print "Files have different bin sizes, exiting..."
    sys.exit()
else :
    print "%i %i %i" %(len(fulllines),len(h1lines),len(h2lines))

fulldata=[]
h1data=[]
h2data=[]
bins=[]

for i in range(len(fulllines)) :
    for n in range(1) : #change to 5 if you want flat bins
        fulldata.append(fulllines[i])
        h1data.append(h1lines[i])
        h2data.append(h2lines[i])
        bins.append(i*5-180+n)

plot(bins,h1data,'b',lw=3,label='1st half')
plot(bins,h2data,'r',lw=3,label='2nd half')
plot(bins,fulldata,'k',lw=3,label='Full simulation')

lg=legend(fontsize=25,labelspacing=0)
lg.get_frame().set_linewidth(0)
title(pdftitle,fontsize=45)
xlim(-180,180)
ylim(0,0.15)
xticks(arange(-180,181,60),fontsize=20)
yticks(fontsize=20)
if grid :
    xlabel(r'$\chi_%s$ dihedral angle' %chi)
    ylabel('Probability')

print "Writing to %s.pdf" %name
savefig("%s.pdf" %name,format='pdf')
