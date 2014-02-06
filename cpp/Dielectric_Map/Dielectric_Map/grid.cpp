//
//  grid.cpp
//  Dielectric_Map
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "grid.h"


void Grid::read_files( std::string & pqrfile, std::string & inputfile, std::string & polarizabilityfile )
{
    pqr.ReadPQR(pqrfile);
    input.ReadInput(inputfile);
    parm.ReadParm(polarizabilityfile);
        
    atoms=pqr.atoms();
    amoeba=parm.amoeba();
    apbsin=input.coord();
}
    
void Grid::grid_diel()
{
    gridindices atom_grid;
    npoints=apbsin.dime[0]*apbsin.dime[1]*apbsin.dime[2];
    
    std::vector<int> nxcounts(npoints), nycounts(npoints),nzcounts(npoints);
    std::vector<double> xshifted_diel(npoints), yshifted_diel(npoints), zshifted_diel(npoints);
    
    for (std::vector<atomdata>::iterator atom=atoms.begin(); atom!=atoms.end();++atom)
    {
        atom_grid=map_atoms_to_grid(*atom);
        map_diel_to_grid(atom_grid.xs, atom_grid.y , atom_grid.z , atom->Alpha, xshifted_diel, nxcounts);
        map_diel_to_grid(atom_grid.x , atom_grid.ys, atom_grid.z , atom->Alpha, yshifted_diel, nycounts);
        map_diel_to_grid(atom_grid.x , atom_grid.y , atom_grid.zs, atom->Alpha, zshifted_diel, nzcounts);
    }
    
    for (int i=0;i<npoints;i++)
    {     
        if (nxcounts[i] > 0){xshifted_diel[i]/=nxcounts[i];}
        if (nycounts[i] > 0){yshifted_diel[i]/=nycounts[i];}
        if (nzcounts[i] > 0){zshifted_diel[i]/=nzcounts[i];}
        if (xshifted_diel[i] == 0){xshifted_diel[i]=apbsin.sdie;}
        if (yshifted_diel[i] == 0){yshifted_diel[i]=apbsin.sdie;}
        if (zshifted_diel[i] == 0){zshifted_diel[i]=apbsin.sdie;}
    }
    xdiel=xshifted_diel;
    ydiel=yshifted_diel;
    zdiel=zshifted_diel;
}

void Grid::map_diel_to_grid(std::vector<int> x, std::vector<int> y, std::vector<int> z, double alpha, std::vector<double> &shifted_diel, std::vector<int> &ncounts )
{
    int grid_point;
    for (std::vector<int>::iterator i=x.begin(); i!=x.end(); ++i )
    {
        for (std::vector<int>::iterator j=y.begin(); j!=y.end(); ++j)
        {
            for (std::vector<int>::iterator k=z.begin(); k!=z.end(); ++k)
            {
                grid_point=indices_to_grid(*i, *j, *k);
                shifted_diel[grid_point] += alpha;
                ncounts[grid_point] += 1;
            }
        }
    }
}

int Grid::indices_to_grid( int x, int y, int z )
{
    if (z >= apbsin.dime[2] || y >= apbsin.dime[1] )
    {
        std::cout << z << " " << y << " " << x << std::endl;
        std::cout << "An error has occurred when mapping atoms to grid points--an atom is off the grid.  This will required more debugging in grid.cpp" << std::endl;
        exit(1);
    }
    return z+y*apbsin.dime[2]+x*apbsin.dime[1]*apbsin.dime[1];
}

// This functions determines which grid points an atom encompasses; it's faster to
// check each atom once than to check each atom dime[0]*dime[1]*dime[2] times.
gridindices Grid::map_atoms_to_grid( atomdata atoms )
{
    gridindices gridpositions;
    double minposition[3];
    double maxposition[3];
    double shiftedmax[3];
    double shiftedmin[3];
    
    // Iterate over x,y,z
    for (int i=0;i<3;i++)
    {
        // Get distance between the atom and the starting point
        minposition[i]=atoms.P[i]-apbsin.min[i]-atoms.R;
        maxposition[i]=atoms.P[i]-apbsin.min[i]+atoms.R;
        shiftedmin[i] =atoms.P[i]-apbsin.min[i]-atoms.R-apbsin.gs[i]/2;
        shiftedmax[i] =atoms.P[i]-apbsin.min[i]+atoms.R-apbsin.gs[i]/2;
        // Divide that distance by the grid spacing to determine how many grid points are in that distance
        minposition[i]/=apbsin.gs[i];
        maxposition[i]/=apbsin.gs[i];
        shiftedmin[i] /=apbsin.gs[i];
        shiftedmax[i] /=apbsin.gs[i];
        // Roud the minimum distance up and the maximum distance down, since the minimum distance is
        // above that grid point and the maximum distance is below that grid point
        minposition[i]= ceil(minposition[i]);
        maxposition[i]=floor(maxposition[i]);
        shiftedmin[i] = ceil(shiftedmin[i]);
        shiftedmax[i] =floor(shiftedmax[i]);
        // If part of the range of grid points is in the box, make sure all of the range is in the box.
        // If it's not, the 'if statement' after this loop will weed those out.
        if (minposition[i] < 0 && maxposition[i] > 0){minposition[i]=0;}
        if (maxposition[i] > apbsin.dime[i]-1 && minposition[i] < apbsin.dime[i]){maxposition[i]=apbsin.dime[i]-1;}
        if (shiftedmin[i]  < 0 && shiftedmax[i] >  0){shiftedmin[i] =0;}
        if (shiftedmax[i]  > apbsin.dime[i]-1 && shiftedmin[i]  < apbsin.dime[i]){shiftedmax[i] =apbsin.dime[i]-1;}
    }
    // Make sure we're only adding in grid points inside the box.
    if (minposition[0] >=0 && minposition[1] >=0 && minposition[2] >=0 &&
        maxposition[0] < apbsin.dime[0] && maxposition[1] < apbsin.dime[1] && maxposition[2] < apbsin.dime[2])
    {
        for (int i=minposition[0];i<=maxposition[0];i++){gridpositions.x.push_back(i);}
        for (int i=minposition[1];i<=maxposition[1];i++){gridpositions.y.push_back(i);}
        for (int i=minposition[2];i<=maxposition[2];i++){gridpositions.z.push_back(i);}
    }
    // And again, make sure we're only including grid points inside the box.
    if (shiftedmin[0] >= 0 && shiftedmin[1] >=0 && shiftedmin[2] >=0 &&
        shiftedmax[0] < apbsin.dime[0] && shiftedmax[1] < apbsin.dime[1] && shiftedmax[2] < apbsin.dime[2])
    {
        for (int i=shiftedmin[0];i<=shiftedmax[0];i++){gridpositions.xs.push_back(i);}
        for (int i=shiftedmin[1];i<=shiftedmax[1];i++){gridpositions.ys.push_back(i);}
        for (int i=shiftedmin[2];i<=shiftedmax[2];i++){gridpositions.zs.push_back(i);}
    }
        /*
         cout << atoms.resid << " " << atoms.name << endl;
         for (int i=0;i<3;i++) {cout<<minposition[i]<< " " << maxposition[i] << endl;}
         cout << endl;
         for (vector<int>::iterator n=gridpositions.xs.begin();n!=gridpositions.xs.end();++n)
         {
         cout << *n << " ";
         }
         cout << endl;
         for (vector<int>::iterator n=gridpositions.ys.begin();n!=gridpositions.ys.end();++n)
         {
         cout << *n << " ";
         }
         cout << endl;
         for (vector<int>::iterator n=gridpositions.zs.begin();n!=gridpositions.zs.end();++n)
         {
         cout << *n << " ";
         }
         cout << endl<<endl;;
         */
    return gridpositions;
}
    
void Grid::write_diel( std::string & outname )
{
    std::string xname, yname, zname;
    dxname(outname,"x.dx",xname);
    dxname(outname,"y.dx",yname);
    dxname(outname,"z.dx",zname);
    write_dx(xname,xdiel,0);
    write_dx(yname,ydiel,1);
    write_dx(zname,zdiel,2);

}

void Grid::dxname( std::string basename, std::string suffix, std::string & dxfilename )
{
    std::string extension=".dx";
    if(basename.rfind(extension) != std::string::npos)
    {
//        std::cout << "Removing extraneous " << dx << " extension from output file name." << std::endl;
        basename.replace(basename.rfind(extension),basename.length(),"");
    }
    dxfilename.append(basename);
    dxfilename.append(suffix);
}

void Grid::write_dx( std::string filename, std::vector<double> grid_diel, int index )
{
    std::vector<double> corner;
    std::string shifted;
    switch ( index )
    {
        case 0 :
            corner.push_back(apbsin.min[0]+apbsin.gs[0]/2);
            corner.push_back(apbsin.min[1]);
            corner.push_back(apbsin.min[2]);
            shifted = "X-SHIFTED";
            break;
        case 1 :
            corner.push_back(apbsin.min[0]);
            corner.push_back(apbsin.min[1]+apbsin.gs[1]/2);
            corner.push_back(apbsin.min[2]);
            shifted="Y-SHIFTED";
            break;
        case 2 :
            corner.push_back(apbsin.min[0]);
            corner.push_back(apbsin.min[1]);
            corner.push_back(apbsin.min[2]+apbsin.gs[2]/2);
            shifted="Z-SHIFTED";
            break;
        default :
            std::cout << "Error writing .dx file.  More debugging likely needed in grid.cpp" << std::endl;
            exit(1);
    }
    
    std::ofstream outputfile;
    std::cout << "Writing to " << filename << std::endl;
    outputfile.open(filename.c_str());
    outputfile << "# Data for APBS 1.3 from Dielectric_Map" << std::endl;
    outputfile << "#"<< std::endl;
    outputfile << "# " << shifted << " DIELECTRIC MAP"<< std::endl;
    outputfile <<"#"<< std::endl;
    outputfile << "object 1 class gridpositions counts " << apbsin.dime[0] << " " << apbsin.dime[1] << " " <<apbsin.dime[2] << std::endl;
    outputfile << std::scientific;
    outputfile << "origin " << corner[0] << " " << corner[1] << " " << corner[2] << std::endl;
    outputfile << "delta " << apbsin.gs[0] << " " << 0.0 << " " << 0.0 << std::endl;
    outputfile << "delta " << 0.0 << " " << apbsin.gs[1] << " " << 0.0 << std::endl;
    outputfile << "delta " << 0.0 << " " << 0.0 << " " << apbsin.gs[2] << std::endl;
    outputfile << "object 2 class gridconnections counts " << apbsin.dime[0] << " " << apbsin.dime[1] << " " <<apbsin.dime[2] << std::endl;
    outputfile << "object 3 class array type double rank 0 items " << npoints << " data follows" << std::endl;
    int col=0;
    for (std::vector<double>::iterator i=grid_diel.begin(); i!=grid_diel.end();++i)
    {
        if (col==2)
        {
            outputfile << *i << std::endl;
            col=0;
        }
        else
        {
            outputfile << *i << " ";
            col++;
        }
    }
    if (col!=0){outputfile << std::endl;}
    outputfile << "attribute \"dep\" string \"positions\"" << std::endl;
    outputfile << "object \"regular positions regular connections\" class field" << std::endl;
    outputfile << "component \"positions\" value 1" << std::endl;
    outputfile << "component \"connections\" value 2" << std::endl;
    outputfile << "component \"data\" value 3" << std::endl;
    outputfile.close();
}