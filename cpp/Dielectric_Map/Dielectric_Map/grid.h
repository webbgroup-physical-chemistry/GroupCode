//
//  grid.h
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#ifndef __Dielectric_Map__grid__
#define __Dielectric_Map__grid__

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "read_pqr.h"
#include "read_input.h"

#endif /* defined(__Dielectric_Map__grid__) */

struct grid
{
    std::vector<double> x, y, z;
};

struct gridindices
{
    std::vector<int> x,y,z;
    std::vector<int> xs,ys,zs;
};



class Grid
{
    // This is for reading files
    CoordinateFile pqr;
    InputFile input;
    ParmFile parm;
    std::vector<atomdata> atoms;
    std::map<std::string,parms> amoeba;
    coords apbsin;
    int npoints;
    
    // This needs to get cleaned up in a little bit
    double x, y, z;
    double grid_point[3], shifted_point[3], xshifted_point[3], yshifted_point[3], zshifted_point[3];
    std::vector<double> xdiel, ydiel, zdiel;
    
public:
    void read_files( std::string & pqrfile, std::string & inputfile, std::string & polarizabilityfile );
    void grid_diel();
    // This functions determines which grid points an atom encompasses; it's faster to
    // check each atom once than to check each atom dime[0]*dime[1]*dime[2] times.
    gridindices map_atoms_to_grid( atomdata atoms );
    // This functions checks to see if a grid point is in any atoms' VDW radii.  This is super slow.
    void map_diel_to_grid( std::vector<int> x, std::vector<int> y, std::vector<int> z, double alpha, std::vector<double> &shifted_diel, std::vector<int> &ncounts );
    int indices_to_grid( int x, int y, int z );
    void write_diel( std::string & outname );
    void dxname( std::string basename, std::string suffix, std::string & dxfilename );
    void write_dx( std::string filename, std::vector<double> grid_diel, int index );

};