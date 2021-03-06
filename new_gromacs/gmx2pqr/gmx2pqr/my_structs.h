//
//  my_structs.h
//  gmx2pqr
//
//  Created by Andrew Ritchie on 5/4/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#ifndef gmx2pqr_my_structs_h
#define gmx2pqr_my_structs_h

#include <stdio.h>
#include <stdlib.h>

#endif

typedef struct gmx2amb gmx2amb;
typedef struct gmx2amb_dat gmx2amb_dat;
typedef struct hbond_array hbond_array;

struct gmx2amb {
    char resname[20];
    char gmxname[20];
    char ambername[20];
    float charge, radius;
};

struct gmx2amb_dat {
    gmx2amb *dat;
    int n;
};

