#! /usr/local/bin/python

import os
import sys
import textwrap

maxwidth = 79

def readFile( filename, kill=True ) :
    try :
        File = open(filename)
        FileLines = File.readlines()
        File.close()
    except :
        print "Cannot open %s"%filename
        if kill :
            print "EXITING..."
            sys.exit()
        return 0
    return FileLines

def printbox( string ) :
    printw( "\n"+"-"*maxwidth )
    printw("|"+" "*(maxwidth-2)+"|")
    for line in textwrap.wrap(string,width=maxwidth-2) :
        printw("|"+line.center(maxwidth-2)+"|")
    printw("|"+" "*(maxwidth-2)+"|")
    printw( "-"*maxwidth +"\n")
    return True

def printcenter( string ) :
    for line in textwrap.wrap(string,width=maxwidth) :
        printw(line.center(maxwidth))
    return True

def printw( string ) :
    lines = string.split('\n')
    for line in lines :
        print textwrap.fill(line,width=maxwidth)
    return True

def backup_outname( filename ) :
    filename = os.path.abspath(filename)
    basename = os.path.basename(filename)
    dirname = os.path.dirname(filename)
    copyname = filename
    n = 1
    if os.path.isfile(filename) :
        while os.path.isfile(copyname) :
            if n == 100 :
                printw( "Will no make 100 copies of %s, exiting..."%filename )
                return
            copyname = "%s/#%s.%i#"%(dirname,basename,n)
            n += 1
        os.rename(filename, copyname)
        printbox( "Backing up %s to %s"%(filename,copyname))
    return filename

def isFile( name, kill=False ) :
    if not os.path.isfile(name) :
        print "\n>>%s not found<<\n"%name
        if kill :
            print "EXITING..."
            sys.exit()
        else : return False
    return True

def writeLines( name, newlines ) :
    if len(name) > 1054 :
        printw("The name is REALLY long!  You probably put the line string there by mistake.  Edit the code if you want this to work...")
        return 0
    backup_outname( name )
    File = open(name,"w")
    File.write(newlines)
    File.close()
    printw("\n>> Done Writing %s <<\n"%name)
    return name

