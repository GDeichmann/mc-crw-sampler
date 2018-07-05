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

#ifndef MC_KERNEL_H
#define MC_KERNEL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void do_mc(t_topology *top,params *p,int natoms,atom *at,
                rvec *x,int **bonds,int nmol,molecule *mol,double ***nonb,double ***nonbf,
                                dihedmap *dhm,int bAppend,char *fn_trjout,int trj_type);

#endif
