//
//  read_pqr.h
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#ifndef __Dielectric_Map__read_pqr__
#define __Dielectric_Map__read_pqr__

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>

#endif /* defined(__Dielectric_Map__read_pqr__) */

struct atomdata
{
    int index, resid;
    std::string mtype, atype, residue, name;
    double Q, R, Q2, Alpha;
    std::vector<double> P;
};

class CoordinateFile
{
    std::string mty, aty, res;
    int ind, rid, Natoms;
    double x, y, z, q, r, q2, alpha;
    std::vector<atomdata> Atoms;
    std::vector<std::string> N;
public:
    std::vector<atomdata> atoms();
    
    int natoms();
    
    std::vector<std::string> n();
    
    void ReadPQR( std::string & filename );
};