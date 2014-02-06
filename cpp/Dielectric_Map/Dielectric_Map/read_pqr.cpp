//
//  read_pqr.cpp
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "read_pqr.h"

using namespace std;

vector<atomdata> CoordinateFile::atoms() {return Atoms;}
int CoordinateFile::natoms() {return Natoms;}
vector<string> CoordinateFile::n() {return N;}

void CoordinateFile::ReadPQR( string & filename )
{
    string line;
    atomdata atom;
    Natoms=0;
    ifstream file (filename.c_str());
    if (file.is_open() )
    {
        while ( file.good() )
        {
            getline( file, line );
            stringstream linestream( line );
            linestream >> mty >> ind >> aty >> res >> rid >> x >> y >> z >> q >> r;
            atom.mtype  = mty;
            atom.index  = ind;
            atom.atype  = aty;
            atom.residue= res;
            atom.resid  = rid;
            atom.P.erase(atom.P.begin(),atom.P.end());
            atom.P.push_back(x);
            atom.P.push_back(y);
            atom.P.push_back(z);
            atom.R      = r;
            atom.Q      = q;
            atom.Q2     = q;
            atom.Alpha  = 1;
            string name;
            name.append(res);
            name.append(".");
            name.append(aty);
            atom.name   = name;
            Atoms.push_back(atom);
            N.push_back(name);
                
                
            Natoms ++;
        }
    }
    else
    {
        std::cout << "Cannot open pqr file " << filename << ". Exiting" << std::endl;
        exit(1);
    }
}
    
