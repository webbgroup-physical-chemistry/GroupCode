//
//  read_input.h
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#ifndef __Dielectric_Map__read_input__
#define __Dielectric_Map__read_input__

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#endif /* defined(__Dielectric_Map__read_input__) */


struct coords
{
    std::vector<int> dime;
    std::vector<double> glen,gcent,min,max,gs;
    double pdie, sdie;
};

struct parameters
{
    std::string atom, residue;
    double q,a;
};

struct parms
{
    double q,a;
};

class InputFile
{
    std::string keyword;
    double value1, value2, value3;
    coords Coord;
    
public:
    coords coord();
    
    void ReadInput( std::string & filename );
};


class ParmFile
{
    std::vector<parameters> parmfile;
    std::string atomname, residuename;
    double charge, polari;
    std::map<std::string, parms> Vars;
    
public:
    std::vector<parameters> parm();
    std::map<std::string, parms> amoeba();
    
    void ReadParm2( std::string & filename );
    // This is like the python ={} object
    void ReadParm( std::string & filename );
};
