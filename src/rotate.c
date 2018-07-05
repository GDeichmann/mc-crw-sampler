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
#include "rotate.h"


void rotate(int n,int *at,rvec axis,double t,rvec origin,rvec *in,rvec *out){ //rodrigues-gibbs formulation 

    int i,atnr;
    rvec src,dest,inc,uax;
    unitv(axis,uax);

    for(i=0;i<n;i++){
        clear_rvec(dest);
        atnr=at[i];
        rvec_sub(in[atnr],origin,src);

        svmul(cos(t),src,inc);
        rvec_inc(dest,inc);

        
        cprod(uax,src,inc);
        svmul(sin(t),inc,inc);
        rvec_inc(dest,inc);

        svmul(iprod(src,uax),uax,inc);
        svmul((1-cos(t)),inc,inc);
        rvec_inc(dest,inc);

        rvec_add(dest,origin,out[atnr]);

    }

    return;

}

void prepare_dih_rot(molecule m,atom *a,int **b,int active,rvec *x,rvec *axis,int *nrota,int **rotlist){ //this will be a rotation around a bond from atom 'Active' to a larger atom

    int i,j=0;
    int nrot=0; //the number of atoms to be modified
    int second;
    int next_in_main;
    int nbonds=0;

    //How many atoms are bonded to active
    for(i=m.first+m.nat;i>active;i--){
        if(b[active][i]==1){
            nbonds++;
            if(nbonds==1){
                next_in_main=i;
            }
        }
    }
   
    int bond;
    int bondcount=0;
    if(nbonds>1){
        bond=rand()%nbonds;
    }
    else{
        bond=0;
    }

    for(i=m.first+m.nat;i>active;i--){
        if(b[active][i]==1){
            j++;
        }
        if(j==bond){
            second=i;
            break;
        }
    }
    
    while(nrot==0 && bondcount < nbonds){ // We do this in case nrot is zero on our first branch

        if(bond==0){ //This is the main branch and therefore easy
            nrot=m.nat+m.first-next_in_main-1;
            second=next_in_main;
        }
        else{ //This is more tricky
            int cycle=0;
            int count=0;
            int start=next_in_main-1;
            while(cycle<bond){ //move down all branches
                start-=count;
                count=0;
                i=start;
                while(i>active){
                    for(j=i;j>=0;j--){
                        if(b[i][j]==1){
                            i=j;
                            count+=1;
                            break;
                        }
                    }
                }
                cycle+=1;
            }
            nrot=count-1;
            second=start-count+1;
        }

        if(nrot==0){
            if(bond==nbonds-1){
                bond=0;
            }
            else{
                bond++;
            }
        }

        bondcount++;
    
    }


    rvec_sub(x[second],x[active],*axis);

    *nrota = nrot;
    *rotlist =(int *) calloc(nrot,sizeof(int));
    for(i=0;i<nrot;i++){
        (*rotlist)[i]=second+i+1;
    }


    return;

}

void prepare_mapped_dih_rot(dihedmap map,molecule m,rvec *x,rvec *axis,int *nrota,int **rotlist,int *act){

    int i;

    int ri1=rand()%map.nactive;
    int active,next,nrot;

    active=map.active[ri1];
    *act=active;

    int ri2=rand()%map.nbond[ri1];

    next=map.next[ri1][ri2];
    nrot=map.nrot[ri1][ri2];

    *nrota = nrot;
    *rotlist =(int *) calloc(nrot,sizeof(int));
    for(i=0;i<nrot;i++){
        (*rotlist)[i]=next+i+1;
    }

    rvec_sub(x[next],x[active],*axis);


    return;
}


