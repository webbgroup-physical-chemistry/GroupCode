//
//  MyOptionParser.cpp
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "MyOptionParser.h"

using namespace std;

vector<string> OptionParser::Options(){return options;}

vector<bool> OptionParser::boolOptions() {return booloptions;}

void OptionParser::assign_header( const string & info )
{
    header=info;
}

void OptionParser::inputparameters()
{
    int i = 0;
    for (vector<string>::iterator n = options.begin(); n != options.end(); ++n)
    {
        cout << "\t" <<  flags[i] << "\n\t\t<" << *n << ">\n\t\t" << helps[i] << endl;
        i++;
    }
    i=0;
    for (vector<bool>::iterator n=booloptions.begin() ; n != booloptions.end(); n++)
    {
        if (*n)
        {
            cout <<"\t" << boolflags[i] << "\n\t\tTrue\n\t\t" << boolhelps[i]<<endl;
        }
        else
        {
            cout <<"\t" << boolflags[i] << "\n\t\tFalse\n\t\t" << boolhelps[i]<<endl;
        }
        i++;
    }
    cout <<"\t" << "-h" << "\n\t\tHelp" << endl;
    return;
};

void OptionParser::userinfo()
{
    std::cout << "\t-p <pqr file>" << std::endl;
    std::cout << "\t-i <input file>" << std::endl;
    std::cout << "\t-f <prm file>" << std::endl;
    std::cout << "\t-o <out file>" << std::endl;
    exit(1);
    return;
};

bool OptionParser::OptionExists(const string & option)
{
    return find(begin,end,option) != end;
};

char * OptionParser::getOption(const string & option)
{
    char ** itr = find(begin,end,option);
    if (itr != end && ++itr !=end)
    {
        return *itr;
    }
    return 0;
};

void OptionParser::read_cmdline(char ** begining, char ** ending)
{
    begin=begining;
    end=ending;
};

void OptionParser::getHelp()
{
    cout << "\n                           Andrew Ritchie                           " << endl;
    cout << "\n" << header << "\n" << endl;
    if ( OptionExists("-h") or OptionExists("--help") )
    {
        inputparameters();
        exit(1);
    }
};

void OptionParser::add_option( const string & option, const string & help, const string & defaultvalue )
{
    flags.push_back(option);
    options.push_back(defaultvalue);
    helps.push_back(help);
};

void OptionParser::add_booloption( const string & option, const string & help, bool defaultvalue )
{
    boolflags.push_back(option);
    booloptions.push_back(defaultvalue);
    boolhelps.push_back(help);
}

void OptionParser::parse_args()
{
    int i=0;
    for (vector<string>::iterator n = flags.begin(); n != flags.end(); ++n)
    {
        if ( OptionExists(*n) )
        {
            options[i]=getOption(*n);
        }
        i++;
    }
    i=0;
    for (vector<string>::iterator n=boolflags.begin(); n != boolflags.end(); ++n)
    {
        if ( OptionExists(*n) )
        {
            if ( booloptions[i] )
            {
                booloptions[i]=false;
            }
            else
            {
                booloptions[i]=true;
            }
        }
        i++;
    }
    getHelp();
    inputparameters();
};
