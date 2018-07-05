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

#include <stdio.h>
#include <stdlib.h>

#include "translate.h"
#include "vec.h"

void translate_mol(molecule m,rvec d,rvec **x){
    
    int i;

    for(i=m.first;i<m.nat+m.first;i++){
        rvec_inc((*x)[i],d);
    }

    return;
}

void get_bead_center(molecule m,rvec *x,atom *a,rvec result){

    int i;
    double w_tot=0,w;
    rvec w_v;
    rvec v_sum;
    clear_rvec(v_sum);

    for(i=m.first;i<m.nat+m.first;i++){
        
        w=a[i].weight;
        svmul(w,x[i],w_v);
        rvec_inc(v_sum,w_v);
        w_tot+=w;

    }

    if(w_tot)
        svmul(1/w_tot,v_sum,result);


    return;
}
