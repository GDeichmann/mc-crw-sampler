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
#include <vec.h>
#include <math.h>

#include "crw_types.h"
#include "read_inp.h"


#include "do_force.h"

double constraint_force(molecule m,rvec c,rvec *f){

    int i;
    rvec total;
    clear_rvec(total);
    double force=0,totw=0;
    rvec e;
    unitv(c,e);

    for(i=m.first;i<m.first+m.nat;i++){
       rvec_inc(total,f[i]);
    }
    force = iprod(e,total);

    return force;

}

void do_nbforce(params *p,rvec *f,int natoms,double ***nb,rvec *x){

    rvec coul,disp,rep,inc;

    int i,j;
    rvec r;
    double inv_dist,inv_sq_dist,inv_6_dist;

    for(i=0;i<natoms;i++){
        clear_rvec(coul);
        clear_rvec(disp);
        clear_rvec(rep);
        for(j=0;j<natoms;j++){
            if(i!=j){
                rvec_sub(x[i],x[j],r);
                inv_sq_dist=1/norm2(r);
		        inv_dist=1/norm(r);
                inv_6_dist=inv_sq_dist*inv_sq_dist*inv_sq_dist;
                unitv(r,r);

                if(p->r_c[0] < 0.0 || inv_sq_dist > 1.0/(p->r_c[0]*p->r_c[0])){
                    svmul(nb[0][i][j]*inv_sq_dist,r,inc);
                    rvec_inc(coul,inc);
                }
                
                if(p->r_c[1] < 0.0 || inv_sq_dist > 1.0/(p->r_c[1]*p->r_c[1])){
                    
                    svmul(-1.0*nb[1][i][j]*inv_6_dist*inv_dist,r,inc);
                    rvec_inc(disp,inc);
                    
                    svmul(nb[2][i][j]*inv_6_dist*inv_6_dist*inv_dist,r,inc);
                    rvec_inc(rep,inc);
                }
                
            }
        }
        copy_rvec(coul,f[i]);
        rvec_inc(f[i],disp);
        rvec_inc(f[i],rep);
    }
    
    return;
}
