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

#ifndef ROTATE_H
#define ROTATE_H

#include <stdio.h>
#include <stdlib.h>
#include <vec.h>
#include <math.h>


void rotate(int n,int *at,rvec axis,double t,rvec origin,rvec *in,rvec *out);

void prepare_dih_rot(molecule m,atom *a,int **bond,int active,rvec *x,rvec *axis,int *nrota,int **rotlist);
void prepare_mapped_dih_rot(dihedmap map,molecule m,rvec *x,rvec *axis,int *nrota,int **rotlist,int *act);

#endif
