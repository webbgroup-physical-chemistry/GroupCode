#! /usr/bin/python

from pylab import *
#from scipy import *
from math import pi, atan2
from matplotlib.pyplot import *
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--table", dest="table", help="Tab delimited data table")
parser.add_option("-f", "--title", dest="title", help="Figure title")
parser.add_option("-o", "--outname", dest="outname", help="<outname>.pdf")
parser.add_option("-a", "--azimuthal", dest="azim", help="If looking at \
the azimuthal angle, enter 'True' to view only -90 to 90 degrees. \
Default=False", default=False)
parser.add_option("-p", "--paper", dest="paper", help="portrait of landscape.  \
Default=portrait", default="portrait")
parser.add_option("-n", "--nitrile", dest="nitrile", help="Looking at a Nitrile? \
Type <-n False> otherwise", default=True)
(options,args) = parser.parse_args()

table = options.table
title = options.title
outname = options.outname
azim = options.azim
paper = str(options.paper)
nitrile = options.nitrile

if outname[-4:] != ".pdf" :
    outname="%s.pdf" %outname

if nitrile != True :
    nitrile=False
                  
def read_table( table ) :
    # Read a tab delimited txt file as a table
    file = open( table, 'r' )
    matrix = {}
    mutants = []
    probes = []
    columns = file.readline().split()
    lines = file.readlines()
    for line in lines :
        row = line.split()
        probe = row[0]
        row = row[1:]
        column = columns[1:]
        if probe not in probes :
            probes.append( probe )
        for i in range( len( row ) ) :
            key = probe + "." + columns[i]
            if columns[i] not in mutants :
                mutants.append( columns[i] )
            matrix[key] = eval( row[i].strip() )
    file.close()
    return matrix, probes, mutants

def grouping_strategy( table, figure_title, outname, azimuthal = False ) :

    matrix, probes, mutants = read_table( table )
    probe_n = len( probes )

    if paper == 'portrait' :
        w = ceil( probe_n/3. )
        h = ceil( probe_n/w )
    else :
        h = ceil( probe_n/3. )
        w = ceil( probe_n/h )
    
    fig = figure(figsize=(8,8))
    suptitle( '%s' %figure_title, fontsize = 20 )
    
    mutant_list = []
    for mutant in mutants :
        if 'e' not in mutant[-1] :
            mutant_list.append(mutant)
    mutant_n = len( mutant_list )


    color = {}
    if 'D30E' in mutant_list : #means Ras and not Rap
        color['WT'] = 'b'
    else :
        color['WT'] = 'k'
    color['D30E_E31K'] = color['D30E+E31K'] = 'k'
    color['D30E'] = color['K31E'] = 'r'
    color['E31K'] = color['E30D'] = 'g'
    color['M'] = 'm'
    color['E30D+K31E'] = color['E30D_K31E'] = 'b'
    color['Monomer'] = 'm'



    index = 1
    for probe in probes :
        c = 0
        ax = fig.add_subplot( h, w, index, polar = True )
        for mutant in mutant_list :
            angles = matrix[probe+"."+mutant]
            error = matrix[probe+"."+mutant+"e"]
            ax.bar( angles*pi/180, (1-(c*.1)), width = ((pi/4)*2*error/45), align='center',\
                    color=color[mutant], alpha = 0.2, fill = True, \
                    edgecolor = color[mutant], ls = 'solid' )
            ax.bar( angles*pi/180, (1-(c*.1)), width = 0.02, align='center', \
                    edgecolor = color[mutant], label = mutant, \
                    color = color[mutant] )
        
            c +=1
        
        index += 1
        xticks( arange(0, .785*360/45, .785*90/45), (r" $0^{o}$","","","") )
        yticks( arange(0), "" )

        if "CNF" in outname :
            plot_title = '$%s_{CN}$' % probe
            print "Looking at cyanophenylalanine... naming appropriately..."
        elif nitrile :
            plot_title = '$%s_{SCN}$' % probe
            print "Looking at thiocyanate... naming appropriately..."
        else :
            plot_title = probe.replace('C','')
            print "Not looking at thiocyanate nor cyanophenylalanine... naming ambiguously..."
        plt.title( '%s' %plot_title, va = 'center' )
        if w <= 2 and h <= 2 :
            subplots_adjust( wspace = -0.4, hspace = .2 )
        if w == 2 and h == 3 :
            subplots_adjust( wspace = -0.4, hspace = .2 )
        if w == 3 and h >= 3 :
            subplots_adjust( wspace = 0, hspace = .2 )
        if w >= 3 and h < 3 :
            subplots_adjust( wspace = 0, hspace = -.4 )
        
    savefig('%s' %outname, format = 'pdf' )

    # Making plot legend
    leg = figure()
    for mutant in mutant_list :
        plot(0,0, label = mutant, color = color[mutant], linewidth = 4 )
        xticks(arange(0),"")
        yticks(arange(0),"")
    legend( loc = 'center' )
    savefig( '%s' %outname.replace('.pdf','_key.pdf'), format = 'pdf' )
    print "%s has been converted to figure %s" %( table, outname )
    return 


def polar2cart( theta ) :
    return cos( theta*pi/180 ), sin( theta*pi/180 )

def cartesian_plot( table, figure_title, outname, azimuthal=False ) :
    matrix, probes, mutants = read_table( table )
    mutant_list = []
    for mutant in mutants :
        if 'e' not in mutant :
            mutant_list.append(mutant)

    color = {}
    color['WT'] = 'r'
    color['D30E'] = color['E30D'] = 'b'
    color['E31K'] = color['K31E'] = 'k'
    color['D30E+E31K'] = color['D30E_E31K'] = \
                         color['E30D+K31E'] = \
                         color['E30D_K31E'] = 'g'

    plotorder = 321
    for probe in probes :
        subplot( plotorder, aspect = 'equal' )
        for mutant in mutant_list :
            angle = matrix[probe+"."+mutant]
            error = matrix[probe+"."+mutant+"e"]
            x, y = polar2cart( angle )
            dx, dy = polar2cart( error )
            plot([0,x],[0,y], color = color[mutant], linewidth = 2)
            

        plot([0,1],[0,0], linewidth = 3, color = 'k' )
        if nitrile :
            print "You indicated that you ARE looking at a nitrile"
            plot_title = '$%s_{SCN}$' % probe
        else :
            plot_title = probe.replace('C','')
            print "You indicated that you are NOT looking at a nitrile"
        plt.title( '%s' %plot_title, va = 'bottom', fontsize = 12 )
        xticks( arange(0), "" )
        yticks( arange(0), "" )

        if azimuthal :
            axis([-1,1,-1,1])
        else:
            axis([-1.05,1.05,-1.05,1.05])

        plotorder += 1

    suptitle( '%s' %figure_title, fontsize = 20 )
    subplots_adjust( wspace = -0.7, hspace = 0.25 )

    # Making plot legend
    leg = figure()
    for mutant in mutant_list :
        plot(0,0, label = mutant, color = color[mutant], linewidth = 4 )
        xticks(arange(0),"")
        yticks(arange(0),"")
    legend( loc = 'center' )
    savefig( '%s_key.pdf' %outname, format = 'pdf' )
    print "%s has been converted to figure %s.pdf" %( table, outname )
    savefig( '%s.pdf' %outname, format = 'pdf' )

    
        
            
#cartesian_plot( 'polartable.txt', 'phtitle', 'phoutname', azimuthal=False )
grouping_strategy( table, title, outname, azimuthal=azim )
#show()

