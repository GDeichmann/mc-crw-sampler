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

#ifndef READ_INP_H
#define READ_INP_H

#include <stdlib.h>
#include <stdio.h>

#include "vec.h"
#include "smalloc.h"
#include "tpxio.h"

void read_atypes(int *n,atype **a);
void read_atypes_top(int *n,atype **a,t_topology *top);
void read_atoms(int *n,atom **a,rvec **x);
void read_atoms_top(int *n,atom **a,t_topology *top,t_inputrec *ir);
void read_coords(char *fn,int n,atom *a,rvec *x);
void read_dihtypes(int *n,dihedtype **dt);
void read_dihed(int *n,dihed **dh);
void read_bonds(int n,int **b);
void read_bonds_top(int **b,t_topology *top);
void read_bl(int n,int **b);
void get_mols(int *n,molecule **m,int natom,atom *a);
void get_mols_top(int *n,molecule **m,t_topology *top);
void init_nonb(params *p,double ***nonb,double ***nonbf,int natoms,atom *at,atype *atp,int **bond);
int nb_count(int a,int b,int **cb);
void premap_dihedrals(molecule m,atom *a,int **b,int **bl,dihedmap *dhm);
int highest_in_branch(int maxat,int active,int **cb);

void init_nonb_from_top(t_topology *top,double ***inter_coeff,double ***nonbf,int natoms,t_inputrec *ir,gmx_mtop_t *mtop,int excl_mode);

#endif
