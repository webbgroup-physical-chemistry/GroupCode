#!/usr/bin/python

from mpmath import exp
from numpy import loadtxt, array, savetxt
import sys

potential=loadtxt( sys.argv[1], "float" )
beta = 0.400907

prob=[]
for potvalue in potential:
    prob.append( exp(-beta*potvalue ) )

prob = array(prob)
prob = prob / prob.sum()

probname = sys.argv[1].replace( "mean", "prob" )
savetxt( probname, prob ) 
