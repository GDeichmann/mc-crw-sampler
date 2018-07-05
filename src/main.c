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
#include <time.h>

#include "vec.h"
#include "tpxio.h"
#include "mtop_util.h"

#include "crw_types.h"
#include "read_inp.h"
#include "do_ener.h"
#include "rotate.h"
#include "mc_kernel.h"
#include "read_par.h"

#include <malloc.h>

int main(int argc,char *argv[]){

    char *fn_conf_in = "conf.gro";
    int append_stat=0;

    char *fn_trjout = "mctraj.gro";
    int trj_type=0; // 0 = .gro; 1 = .xtc

    char *fn_tprin = "topol.tpr";
    int excl_mode = 0;
    
    parse_cmdline(argc,argv,&fn_conf_in,&fn_trjout,&trj_type,&append_stat,&fn_tprin,&excl_mode);

    t_inputrec *ir = calloc(1,sizeof(t_inputrec));
    t_state *state = calloc(1,sizeof(t_state));
    gmx_mtop_t *mtop = calloc(1,sizeof(gmx_mtop_t));
    
    read_tpx_state(fn_tprin,ir,state,NULL,mtop);
    t_topology top = gmx_mtop_t_to_t_topology(mtop);

    srand(time(NULL));

    int i,j,k;

    table tab[3];
   /* 
    for(i=0;i<3;i++){
        tab[i].nrow=100000;
        tab[i].prec=10000;
        tab[i].val=(double *) calloc(100000,sizeof(double));
        for(j=1;j<100001;j++){
            if(i==0)
                tab[i].val[j-1]=1/( (((double)j)-0.5) /10000);
            else
                tab[i].val[j-1]=1/pow(( (((double)j)-0.5) /10000),6*i);
        }
    }*/


    params p;
    
    fprintf(stderr,"Reading parameters...\n");
    read_param_ir(&p,ir);


    int natoms;
    atom *at;
    rvec *x;
    
    fprintf(stderr,"Reading atoms...\n");
    read_atoms_top(&natoms,&at,&top,ir);
    x = (rvec *) calloc(natoms,sizeof(rvec));
    read_coords(fn_conf_in,natoms,at,x);

    if(natoms != top.atoms.nr){
        fprintf(stderr,"FATAL ERROR: Number of atoms in %s not equal number of atoms in data.atoms\n",fn_tprin);
        exit(2);
    }

    int **bonds;

    bonds = (int **) calloc(natoms,sizeof(int *));
    
    for(i=0;i<natoms;i++){
        bonds[i] = (int *) calloc(natoms,sizeof(int));
    }

    fprintf(stderr,"Reading bonds...\n");
    read_bonds_top(bonds,&top);

    int nmol;
    molecule *mol;

    fprintf(stderr,"Reading molecules...\n");
    get_mols_top(&nmol,&mol,&top);

    double ***nonb;
    double ***nonbf;

    nonb = (double ***) calloc(3,sizeof(double **)); //0: coulomb; 1: LJ-6; 2: LJ-12
    nonbf = (double ***) calloc(3,sizeof(double **));
 
    for(j=0;j<3;j++){
        nonb[j] = (double **) calloc(natoms,sizeof(double *));
        nonbf[j] = (double **) calloc(natoms,sizeof(double *));
    
        for(i=0;i<natoms;i++){
            nonb[j][i] = (double *) calloc(natoms,sizeof(double));
            nonbf[j][i] = (double *) calloc(natoms,sizeof(double));
        }
    }

    fprintf(stderr,"Reading nonb parameters...\n");
    init_nonb_from_top(&top,nonb,nonbf,natoms,ir,mtop,excl_mode);

    dihedmap *dihmap;
    dihmap = (dihedmap *) calloc(nmol,sizeof(dihedmap));
 
    int **blacklist; //Works like the bond matrix and will prevent certain bonds from becoming active for dihedal rotation
    blacklist = (int **) calloc(natoms,sizeof(int *));
    
    for(i=0;i<natoms;i++){
        blacklist[i] = (int *) calloc(natoms,sizeof(int));
    }
    fprintf(stderr,"Premapping dihedral rotations ...\n");
    read_bl(natoms,blacklist);


    for(i=0;i<nmol;i++){
        premap_dihedrals(mol[i],at,bonds,blacklist,&(dihmap[i]));

        for(j=0;j<dihmap[i].nactive;j++){
            printf("Active %i Nbonds: %i \n",dihmap[i].active[j]+1,dihmap[i].nbond[j]);
            for(k=0;k<dihmap[i].nbond[j];k++){
                printf("\t%i %i\n",dihmap[i].next[j][k]+1,dihmap[i].nrot[j][k]);
            }

        }
    }
    

    for(i=0;i<natoms;i++){
        printf("%f %f %f\n",x[i][0],x[i][1],x[i][2]);
    }


    do_mc(&top,&p,natoms,at,x,bonds,nmol,mol,nonb,nonbf,dihmap,append_stat,fn_trjout,trj_type);


    free(at);
    free(x);
    free(mol);
    free(nonb);
    free(nonbf);
    free(dihmap);
    free(blacklist);

    fprintf(stderr,"\nRien ne va plus!\n\n");


    return 0;
}
