#include "gmx_tilt.h"

/* Proposed workflow : 
    1. Remove solvent from .gro
    2. Remove solvent from .xtc
    3. VMD Stamp Structural Alignment of GTPase, NoSol.gro to 1LFD Chain B
    4. trjconv align trajectory to aligned gro
    5. Make 1LFD into a .gro file
    6. Make index file of Ral for alignment
*/

int main(int argc, char * argv[])
{
    std::cout << "Hello World!" << std::endl;
    gmx_tilt(argc, argv);
    return 0;
}
