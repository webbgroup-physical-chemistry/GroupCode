//
//  read_input.cpp
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "read_input.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

coords InputFile::coord() {return Coord;}
    
void InputFile::ReadInput( string & filename )
{
    string line;
    ifstream file (filename.c_str());
    if (file.is_open() )
    {
        while ( file.good() )
        {
            getline( file, line );
            stringstream linestream( line );
            linestream >> keyword >> value1 >> value2 >> value3;
            if ( keyword.compare("dime") == 0 )
            {
                Coord.dime.push_back(value1);
                Coord.dime.push_back(value2);
                Coord.dime.push_back(value3);
            }
            if (keyword.compare("gcent") == 0 )
            {
                Coord.gcent.push_back(value1);
                Coord.gcent.push_back(value2);
                Coord.gcent.push_back(value3);
            }
            if (keyword.compare("glen") == 0 )
            {
                Coord.glen.push_back(value1);
                Coord.glen.push_back(value2);
                Coord.glen.push_back(value3);
            }
            if (keyword.compare("sdie") == 0 )
            {
                Coord.sdie=value1;
            }
            if (keyword.compare("pdie") == 0 )
            {
                Coord.pdie=value2;
            }
        }
    }
    else
    {
        std::cout << "Cannot open apbs input file " << filename << ". Exiting" << std::endl;
        exit(1);
    }
    
    // boxes are going to be shifted by + glenC/dimeC/2
    for (int i=0;i<3;i++)
    {
        Coord.min.push_back(Coord.gcent[i]-Coord.glen[i]/2);
        Coord.max.push_back(Coord.gcent[i]+Coord.glen[i]/2);
        Coord.gs.push_back(Coord.glen[i]/(Coord.dime[i]-1));
    }
    
//    cout << "Bottom corner: \t( " << Coord.min[0] << " " << Coord.min[1] << " " << Coord.min[2] << " )" << endl;
//    cout << "Top corner: \t (" << Coord.max[0] << " " << Coord.max[1] << " " << Coord.max[2] << " )" << endl;
}


vector<parameters> ParmFile::parm() {return parmfile;}
map<string, parms> ParmFile::amoeba() {return Vars;}
    
void ParmFile::ReadParm2( string & filename )
{
    parameters parmdata;
    string line;
    ifstream file (filename.c_str());
    if (file.is_open() )
    {
        while ( file.good() )
        {
            getline( file, line );
            stringstream linestream( line );
            linestream >> residuename >> atomname >> charge >> polari;
            parmdata.residue=residuename;
            parmdata.atom=atomname;
            parmdata.q=charge;
            parmdata.a=polari;
            parmfile.push_back(parmdata);
        }
    }
}
// This is like the python ={} object
void ParmFile::ReadParm( string & filename )
{
    string line;
    ifstream file (filename.c_str());
    if (file.is_open() )
    {
        while ( file.good() )
        {
            getline( file, line );
            stringstream linestream( line );
            linestream >> residuename >> atomname >> charge >> polari;
            string str;
            str.append(residuename);
            str.append(".");
            str.append(atomname);
            Vars[str].q=charge;
            Vars[str].a=polari;
        }
    }
    else
    {
        std::cout << "Cannot open force field parameter file " << filename << ". Exiting" << std::endl;
        exit(1);
    }
}

