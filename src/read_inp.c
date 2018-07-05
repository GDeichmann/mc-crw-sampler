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

#include <stdlib.h>
#include <stdio.h>

#include "vec.h"
#include "smalloc.h"
#include "crw_types.h"
#include "read_inp.h"

void read_atoms_top(int *n,atom **a,t_topology *top,t_inputrec *ir){
//read atom poperties and pull weights from top and ir

    int i,j;

    *n = top->atoms.nr;
    *a = (atom *) calloc(*n,sizeof(atom));

    for(i=0;i<*n;i++){

        (*a)[i].m = top->atoms.atom[i].m;
        (*a)[i].q = top->atoms.atom[i].q;
        (*a)[i].type = top->atoms.atom[i].type;

        (*a)[i].name=(char *) calloc(6,sizeof(char));
        (*a)[i].rname=(char *) calloc(6,sizeof(char));

    }

    t_pull *p= ir->pull;
//First check whether we can do CRW here
    if(p == NULL || ir->ePull != epullCONSTRAINT){
        fprintf(stderr,"FATAL ERROR: No com-pulling section in your .tpr file or pull-mode is not set to 'constraint' %i \n",p->eGeom);
        exit(4);
    }

    if(p->ngrp != 1){
        fprintf(stderr,"FATAL ERROR: There are %i pull groups. There should be exactly one!\n",p->ngrp);
        exit(4);
    }
//Then loop over groups and get weight for COM
    double pull_weight;
    for(i=0;i<p->ngrp+1;i++){
//        printf("Group %i has %i atoms\n",i,p->grp[i].nat);
        for(j=0;j<p->grp[i].nat;j++){
            if(p->grp[i].nweight){
                pull_weight=p->grp[i].weight[j];
            }
            else{
                pull_weight=1.0;
            }

            int pull_atom= p->grp[i].ind[j];
            //Now assign to atom structure or whatever
            (*a)[pull_atom].weight = (*a)[pull_atom].m*pull_weight;
        }
    }

    return;
}

void read_coords(char *fn,int n,atom *a,rvec *x){

    int i,j;
    size_t length=200;
    char *line = (char *) calloc(length,sizeof(char));
    char resnum[6];
    resnum[5] = '\0';

    FILE *in = fopen(fn,"r");

    if(!in){
        fprintf(stderr,"FATAL ERROR: Can not open conf. file %s\n",fn);
        exit(1);
    }

    getline(&line,&length,in); //one empty line as comment
    
    getline(&line,&length,in);
    sscanf(line,"%i",&i);

    if(i!=n){
        fprintf(stderr,"FATAL ERROR: %i atoms in conf. file and %i atoms in atom definition\n",i,n);
        exit(1);
    }

    for(i=0;i<n;i++){
        getline(&line,&length,in);

        for(j=0;j<5;j++){
            resnum[j] = line[j];
            a[i].rname[j] = line[j+5];
            a[i].name[j] = line[j+10];
        }
        
        a[i].rnum = (int) atol(resnum);

        sscanf(&(line[20]),"%f %f %f",&(x[i][0]),&(x[i][1]),&(x[i][2]));

    }


    fclose(in);
    free(line);

    return;
}

void read_bonds_top(int **b,t_topology *top){

    int i,j,nbonds=0;
    t_idef id=top->idef;

    int at1,at2;
    for(i=F_BONDS;i<=F_TABBONDS;i++){
        for(j=0;j<id.il[i].nr;){

            at1=((( id.il[i].iatoms[j+1]  )<( id.il[i].iatoms[j+2]  ))?( id.il[i].iatoms[j+1]  ):( id.il[i].iatoms[j+2] ));
            at2=((( id.il[i].iatoms[j+1]  )>( id.il[i].iatoms[j+2]  ))?( id.il[i].iatoms[j+1]  ):( id.il[i].iatoms[j+2] ));
            b[at1][at2] = 1;
            b[at2][at1] = 1;
            j+=3;
            nbonds++;
        }
    }

    return;
}

void read_bl(int n,int **b){ 

    int i,j,nbonds;
    int at1,at2;
    size_t length=200;
    char *line = (char *) calloc(length,sizeof(char));

    FILE *in = fopen("data.blacklist","r");
    if(!in){
        free(line);
        return;
    }
    getline(&line,&length,in); //one empty line as comment
    
    getline(&line,&length,in);
    sscanf(line,"%i",&nbonds);

    for(i=0;i<nbonds;i++){
        getline(&line,&length,in);
        sscanf(line,"%i %i %i",&j,&at1,&at2); 
        b[at1-1][at2-1]=1;
        b[at2-1][at1-1]=1;
    }

    free(line);
    fclose(in);
    return;
}

void get_mols_top(int *n,molecule **m,t_topology *top){

    int i;
    *n = top->mols.nr;
    *m=(molecule *) calloc(*n,sizeof(molecule));
    for(i=1;i<*n+1;i++){
        (*m)[i-1].first = top->mols.index[i-1];
        (*m)[i-1].nat = top->mols.index[i] - (*m)[i-1].first;
    }

    return;

}


void init_nonb_from_top(t_topology *top,double ***inter_coeff,double ***nonbf,int natoms,t_inputrec *ir,gmx_mtop_t *mtop,int excl_mode){
//This function prepares the interaction matrix for the nonbonded interactions. And works in 4 steps
//1. Get the regular LJ and Coulomb parameters based on LJ type and charge
//2. Exclude according to exclusions from .tpr file (all coeffs are set to zero)
//3. Get the modified LJ-14 and Q-14 interactions, if present !!IMPORTANT: Coulomb-14 will only be set if there are LJ-14 interactions
//to make sure this is the case you may check your .tpr file with gmxdump
//4. Exclude based on the energygrp-exclusion in the parameter file

    int i,j;
    const double es_const = 138.935485;
    int ljtype;
    double q1,q2;
//Part 1.
    for(i=0;i<natoms;i++){
        for(j=i+1;j<natoms;j++){
            
            ljtype = top->atoms.atom[i].type*top->atomtypes.nr + top->atoms.atom[j].type;
            q1=top->atoms.atom[i].q;
            q2=top->atoms.atom[j].q;
            
            inter_coeff[0][i][j] = q1*q2*es_const;
            inter_coeff[1][i][j] = top->idef.iparams[ljtype].lj.c6;
            inter_coeff[2][i][j] = top->idef.iparams[ljtype].lj.c12;
            
            inter_coeff[0][j][i] = inter_coeff[0][i][j];
            inter_coeff[1][j][i] = inter_coeff[1][i][j];
            inter_coeff[2][j][i] = inter_coeff[2][i][j];
        }                                                                                                                            
    }

//Part 2.
    int k;
    for(i=0;i<top->excls.nr;i++){
        for(j=top->excls.index[i];j<top->excls.index[i+1];j++){
            k=top->excls.a[j];
            
            inter_coeff[0][k][i] = 0.0;
            inter_coeff[1][k][i] = 0.0;
            inter_coeff[2][k][i] = 0.0;

            inter_coeff[0][i][k] = 0.0;
            inter_coeff[1][i][k] = 0.0;
            inter_coeff[2][i][k] = 0.0;
        }
    }

//Part 3.
    i=0;
    int type14;
    double fudgeqq = top->idef.fudgeQQ;
    
    while(i<top->idef.il[F_LJ14].nr){
        
        type14 = top->idef.il[F_LJ14].iatoms[i];
        int at1 = top->idef.il[F_LJ14].iatoms[i+1];
        int at2 = top->idef.il[F_LJ14].iatoms[i+2];
        
        
        inter_coeff[0][at1][at2] = top->atoms.atom[at1].q*top->atoms.atom[at2].q*es_const*fudgeqq;
        inter_coeff[1][at1][at2] = top->idef.iparams[type14].lj14.c6A;
        inter_coeff[2][at1][at2] = top->idef.iparams[type14].lj14.c12A;
        
        inter_coeff[0][at2][at1] = top->atoms.atom[at1].q*top->atoms.atom[at2].q*es_const*fudgeqq;
        inter_coeff[1][at2][at1] = top->idef.iparams[type14].lj14.c6A;
        inter_coeff[2][at2][at1] = top->idef.iparams[type14].lj14.c12A;
        
        i+=3;
    }


//Part 4
    //check whether there are energy groups
    int ngener = ir->opts.ngener;
    if(ngener>1){ //There will always be 1 "rest" energy group 
        printf("There are %i energy groups\n",ngener);
    }
    //cross-check with topology
    if(ngener != mtop->groups.grps[egcENER].nr){
        printf("Something's fucky %i %i\n",ngener,mtop->groups.grps[egcENER].nr);
    }

    //check whether there are exclusions
    int nexcl_pair=0;
    for(i=0;i<ngener;i++){
        for(j=i;j<ngener;j++){
            nexcl_pair += (ir->opts.egp_flags[i*ngener+j] & EGP_EXCL);
        }
    }
    int l; 
    //collect atoms in exclusion groups
    printf("There are %i exclusion pairs\n",nexcl_pair);
    for(i=0;i<ngener;i++){
        for(j=i;j<ngener;j++){
            if (ir->opts.egp_flags[i*ngener+j] & EGP_EXCL){
                int grp1 = i;
                int grp2 = j;
                for (k=0;k<natoms;k++){
                    for (l=0;l<natoms;l++){
                        if(mtop->groups.grpnr[egcENER][k] == grp1  && 
                            mtop->groups.grpnr[egcENER][l] == grp2){
                           
                            if(excl_mode == 0 || excl_mode == 1){
                                inter_coeff[0][k][l] = 0.0;
                                inter_coeff[0][l][k] = 0.0;
                            }

                            if(excl_mode == 0 || excl_mode == 2){
                                inter_coeff[1][k][l] = 0.0;
                                inter_coeff[2][k][l] = 0.0;
                                inter_coeff[1][l][k] = 0.0;
                                inter_coeff[2][l][k] = 0.0;
                            }
                            

                        }
                    }
                }
            }
        }            
    }
    
    for(i=0;i<natoms;i++){
        for(j=0;j<natoms;j++){
            nonbf[0][i][j] = inter_coeff[0][i][j];
            nonbf[1][i][j] = inter_coeff[1][i][j]*6.0;
            nonbf[2][i][j] = inter_coeff[2][i][j]*12.0;
        }
    }

    return;
}

int nb_count(int a,int b,int **cb){

    int count=0,current;
    int swap=0,i;

    if(a==b){
        return 0;
    }

    if(a<b){
        swap=b;
        b=a;
        a=swap;        
    }
    
    current=a;

    while(current>b){ //We move down the main chain starting at a
        for(i=current;i>=0;i--){
            if(cb[current][i]==1){
                current=i;
                count+=1;
                break;
            }
        }
    }
//Now we are at the branching point where the branch of b hits the main chain

    int branchpoint=current;
    current=b;

    while(current>branchpoint){ //If b is in a branch we move down this branch
        for(i=current;i>=0;i--){
            if(cb[current][i]==1){
                current=i;
                count+=1;
                break;
            }
        }
    }

    return count;

}

void premap_dihedrals(molecule m,atom *a,int **b,int **bl,dihedmap *dhm){

    int i,j,k,active;
    int **next;
    int **nrot;
    int *nbond;
    int nbonds;
    int atom;
    nbond =(int*) calloc(m.nat,sizeof(int));
    next = (int**) calloc(m.nat,sizeof(int*));
    nrot = (int**) calloc(m.nat,sizeof(int*));

    for(active=0;active<m.nat;active++){
        atom=m.first+active;

        nbonds=0;

        if(atom == m.first+m.nat){
            nbond[active]=0;
            break;
        }

        for(i=m.first+m.nat;i>atom;i--){
            if(b[atom][i]==1){
                nbonds++;
            }
        }

        nbond[active] = nbonds;
        
        if(nbonds){
            next[active] = (int*) calloc(nbonds,sizeof(int));
            nrot[active] = (int*) calloc(nbonds,sizeof(int));
        }

        nbonds=0;
        for(i=m.first+m.nat;i>atom;i--){
            
            if(b[atom][i]==1){
                nbonds++;
                next[active][nbonds-1]=i;
            }
        }
        
        for(k=0;k<nbond[active];k++){

            if(bl[atom][next[active][k]]){//First check whether this bond is blacklisted
                nrot[active][k]=0;
                break;
            }
            
            if(k==0){ //We rotate all up to the highest in our respective branch
                nrot[active][k]=highest_in_branch(m.first+m.nat-1,active+m.first,b)-next[active][k];
            }
            else{ //This is more tricky
                int cycle=0;
                int count=0;
                int start=next[active][k-1]-1;

                nrot[active][k]=next[active][k-1]-next[active][k]-1;
                
            }

//            printf("\tBond:%i Next:%i Nrot: %i\n",k+1,next[active][k]+1,nrot[active][k]);
        }

    } 

    int non_zero=0;
    int non_zero_bonds[m.nat];

    for(active=0;active<m.nat;active++){ //Do some counting here
        non_zero_bonds[active]=0;
        if(nbond[active]){
            non_zero++;
            for(i=0;i<nbond[active];i++){
                if(nrot[active][i]){
                    non_zero_bonds[active]++;
                }
            }
            if(!non_zero_bonds[active]){
                non_zero--;
            }
        }
    }
   
    dhm->nactive=non_zero;
    dhm->active = (int*) calloc(non_zero,sizeof(int));
    dhm->nbond = (int*) calloc(non_zero,sizeof(int));
    dhm->next = (int**) calloc(non_zero,sizeof(int*));
    dhm->nrot = (int**) calloc(non_zero,sizeof(int*));

    i=0;

    for(active=0;active<m.nat;active++){ //Actual value assignment is done here
        if(non_zero_bonds[active]){
            dhm->active[i]=active+m.first;
            dhm->nbond[i]=non_zero_bonds[active];
            dhm->next[i]=(int*) calloc(non_zero_bonds[active],sizeof(int));
            dhm->nrot[i]=(int*) calloc(non_zero_bonds[active],sizeof(int));

            j=0;
            for(k=0;k<nbond[active];k++){
                if(nrot[active][k]){
                    dhm->next[i][j]=next[active][k];
                    dhm->nrot[i][j]=nrot[active][k];
                    j++;
                }
            }

            i++;
        }
    }

    free(nbond);
    free(next);
    free(nrot);
    
    return;
}

int highest_in_branch(int maxat,int active,int **cb){

    int count=0,current;
    int i,j;
    
    for(j=maxat;j>active;j--){

        current=j;

        while(current>active){ //We move down the main chain starting at a
            for(i=current;i>=0;i--){
                if(cb[current][i]==1){
                    current=i;
                    count+=1;
                    break;
                }
            }
        }
    //Now we are at the branching point where the branch of b hits the main chain

        int branchpoint=current;

        if(branchpoint==active){
            return j;
        }
        
    }

    return active;

}

