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

#ifndef DO_ENER_H
#define DO_ENER_H

#include <stdio.h>
#include <stdlib.h>

#include "vec.h"
#include "maths.h"
#include "types/idef.h"
//#include "crw_types.h"

double do_bonded(int min,int max,rvec *x,t_idef *id);
double dihed_pot(dihed dh,dihedtype *type,rvec *x);

double mol_intra_pot(molecule m,dihed *dh,dihedtype *dhtp,rvec *x);

double do_nonb(params *p,int natoms,double ***nb,rvec *x);
double do_nonb_tab(table *t,params *p,int natoms,double ***nb,rvec *x);

#endif
