//
//  main.cpp
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "MyOptionParser.h"
#include "grid.h"


using namespace std;

int main(int argc, char * argv[])
{
    /* set defaults */
    Grid* dx = new Grid;
    
    OptionParser parser;
    parser.assign_header("This is intended to build a dielectric map for use in APBS.  It uses AMOEBA\n\
                         point charges and polarizabilities from a TINKER .prm file, an existing .pqr\n\
                         file, and an existing .in file");
    parser.read_cmdline(argv, argv+argc);
    //parser.add_option("-p", "Name of the coordinate .pqr file", "/Users/ritchie/Desktop/tmp_apbs/dielmap/7-WT+G28C_180-180-frame80.pqr");
    parser.add_option("-p", "Name of the coordinate .pqr file", "frame.pqr");
    //parser.add_option("-i", "Name of the APBS input .in file", "/Users/ritchie/Desktop/tmp_apbs/dielmap/testdiel.in");
    parser.add_option("-i", "Name of the APBS input .in file", "apbs.in");
    //parser.add_option("-f", "Name of the charge and polarization mapping file.", "/Users/ritchie/Desktop/tmp_apbs/dielmap/Polarizability_mapping.DAT");
    parser.add_option("-f", "Name of the charge and polarization mapping file.", "FF.DAT");
    parser.add_option("-o", "Dielectric map output file name", "diel.dx");
    parser.add_booloption("-v", "Be Verbose", false);
    parser.parse_args();
    string pqrfile  = parser.Options()[0];
    string infile   = parser.Options()[1];
    string polfile  = parser.Options()[2];
    string outfile  = parser.Options()[3];
    bool bVerbose   = parser.boolOptions()[0];
    
    if (bVerbose){std::cout << "The verbose flag currently does nothing" << std::endl;}
    
    
    dx->read_files(pqrfile,infile,polfile);
    dx->grid_diel();
    dx->write_diel(outfile);
    
    return 0;
}

