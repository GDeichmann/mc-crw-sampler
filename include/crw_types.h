/*
  Copyright (C) 2018   Gregor Deichmann (deichmann@cpc.tu-darmstadt.de)
  
  This file is part of mc-crw-sampler.
  
  mc-crw-sampler is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
  
  mc-crw-sampler is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef CRW_TYPES_H
#define CRW_TYPES_H

#include <stdio.h>
#include <stdlib.h>


#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923
#endif

#ifndef M_2PI
#define M_2PI       6.28318530717958647692
#endif

#ifndef M_SQRT2
#define M_SQRT2 sqrt(2.0)
#endif

#ifndef M_1_PI
#define M_1_PI      0.31830988618379067154
#endif

typedef struct {

    int nactive;
    int *active;
    int *nbond; //bonds per active atom
    int **next; 
    int **nrot;

} dihedmap;

typedef struct {

    int nrow;
    int prec;
    double *val;

} table;

typedef struct {

    long unsigned int nstep; //uint 1
    double kT; //! ref-t
    double beta; //! idem
    int trajfreq; //uint2
    int logfreq; //uint3
    int flogfreq; //uint4
    double p_rot;  //ureal1
    int crw_type; //! 
    int comb_rule; //! same as in gromacs: 1: C6-C12  2: Lorentz Berthelot (sig+eps); 3: OPLS-type all geometric (sig+eps)
    double r_c[2]; //!  0 is coulomb; 1 is vdw
    int nrexcl; //!
    double max_rot[2]; //ureal2+3 index equals move nr
    int restart;//! 1: try restart run

} params;

typedef struct {

    float m;
    float q;
    double sig;
    double eps;
    
} atype;

typedef struct {

    int type;
    float weight;
    float m;
    float q;
    char *name;
    int rnum;
    char *rname;

} atom;


typedef struct {

    int at1,at2,at3,at4;
    int type;

} dihed;

typedef struct {
    
    double params[6];

} dihedtype;

typedef struct {

    int nat;
    int first;

} molecule;

#endif

