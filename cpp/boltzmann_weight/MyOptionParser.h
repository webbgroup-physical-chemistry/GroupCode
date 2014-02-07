//
//  MyOptionParser.h
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#ifndef __Dielectric_Map__MyOptionParser__
#define __Dielectric_Map__MyOptionParser__

#include <algorithm>
#include <vector>
#include <iostream>
#include <stdio.h>

#endif /* defined(__Dielectric_Map__MyOptionParser__) */

class OptionParser
{
private:
    char ** begin, ** end;
    std::string header;
    std::vector<std::string> flags, options, helps, boolflags, boolhelps;
    std::vector<bool> booloptions;
    
public:
    std::vector<std::string> Options();
    std::vector<bool> boolOptions();
    
    void assign_header( const std::string & info );
    
    void inputparameters();
    
    void userinfo();
    
    bool OptionExists(const std::string & option);
    
    char * getOption(const std::string & option);
    
    void read_cmdline(char ** begining, char ** ending);
    
    void getHelp();
    
    void add_option( const std::string & option, const std::string & help, const std::string & defaultvalue );
    
    void add_booloption( const std::string & option, const std::string & help, bool defaultvalue );
    
    void parse_args();
};



