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
#include "xtcio.h"
#include "gmxfio.h"

#include "crw_types.h"
#include "do_ener.h"
#include "do_force.h"
#include "rotate.h"
#include "translate.h"
#include "mc_kernel.h"

#define BOX_SIZE 10.0
#define XTC_PREC 1000

void get_continuation(int trj_type,const char *fn_mclog,const char *fn_trj,const char *fn_flog,rvec *x_out,int natoms,unsigned long int *start_sim){
//looks up the last number in files mclog and flog and compares for consistency
//for .gro trajectories this will probably not work (not necessarily step information present) so we compare natoms and get the last frame as x 

    unsigned long int end_fl,end_ml;

    FILE *ml = fopen(fn_mclog,"r");
    if(ml == NULL){
        fprintf(stderr,"Could not open file %s for appending\n",fn_mclog);
        exit(3);
    }

    FILE *tr = NULL;
    t_fileio *trj = NULL;
    if(trj_type == 0){
        tr = fopen(fn_trj,"r");
        if(tr == NULL){
            fprintf(stderr,"Could not open file %s for appending\n",fn_trj);
            exit(3);
        }
    }
    else if(trj_type == 1){
        trj = open_xtc(fn_trj,"r");
    }

    FILE *fl = fopen(fn_flog,"r");
    if(fl == NULL){
        fprintf(stderr,"Could not open file %s for appending\n",fn_flog);
        exit(3);
    }
    
    size_t linelen=1024;
    char *line = (char *) calloc(linelen,sizeof(char));

    while(!feof(fl)){
        getline(&line,&linelen,fl);
    }
    end_fl = atol(line);

    while(!feof(ml)){
        getline(&line,&linelen,ml);
    }
    end_ml = atol(line);

    if(end_ml != end_fl){
        fprintf(stderr,"FATAL ERROR: Last entry to %s is at %li steps, while last entry to %s is at %li steps\n",fn_mclog,end_ml,fn_flog,end_fl);
        exit(3);
    }

    *start_sim = end_fl;

    unsigned long int nframe,nline=0;
    long int frame_start;

    if(x_out != NULL){
        if(trj_type == 0){
            while(!feof(tr)){
                if(nline%(natoms+3) == 1  && !feof(tr)){ //this gives us the second line but that's fine 
                    frame_start = ftell(tr);
                }
                getline(&line,&linelen,tr);
                nline++;
                if(nline == 2){
                    if(atoi(line) != natoms){
                        fprintf(stderr,"Number of atoms in trajectory %i, while Number of atoms in input %i!\n",atoi(line),natoms);
                    }
                }
            }
            nframe = nline / (natoms+3);

            fseek(tr,frame_start,SEEK_SET);
            
            getline(&line,&linelen,tr);
            int i;
            for(i=0;i<natoms;i++){
                getline(&line,&linelen,tr);

                sscanf(&(line[20]),"%f %f %f",&(x_out[i][0]),&(x_out[i][1]),&(x_out[i][2]));
            }
        }
        else if(trj_type == 1){
            matrix box;
            real prec;
            int step=0;
            real time;
            gmx_bool bOk;
            //read_first_xtc(trj,&natoms,&step,&time,box,&x_out,&prec,&bOk);
            xtc_seek_frame(trj,end_fl,natoms);

            read_next_xtc(trj,natoms,&step,&time,box,x_out,&prec,&bOk);
        }
    }

    free(line);
    fclose(ml);
    fclose(fl);
    if(trj_type == 0 && tr != NULL){
        fclose(tr);
    }
    if(trj_type == 1 && trj != NULL){
        close_xtc(trj);
    }
    return;
}

void init_mc(t_idef *id,params *p,double *bond_pot,double *tot_bond,double *nonb_pot,rvec *ini_pos,int natoms,atom *at,                
                rvec *x,int nmol,molecule *mol,double ***nonb){

    int i;
    *tot_bond =0;

    for(i=0;i<nmol;i++){
        bond_pot[i]=do_bonded(mol[i].first,mol[i].first+mol[i].nat,x,id);
        *tot_bond+=bond_pot[i];
        get_bead_center(mol[i],x,at,ini_pos[i]);
    }
    
    *nonb_pot = do_nonb(p,natoms,nonb,x);
    

    return;
}

void write_gro_frame(FILE *traj,unsigned long int step,int natoms,atom *at,rvec *x){

    fprintf(traj,"So long and thanks for all the fish; step= %li\n",step);
    fprintf(traj,"%i\n",natoms);

    int i;
    for(i=0;i<natoms;i++){
        fprintf(traj,"%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n",at[i].rnum,at[i].rname,at[i].name,i+1,x[i][0],x[i][1],x[i][2]);
    }

    fprintf(traj,"%g %g %g\n",BOX_SIZE,BOX_SIZE,BOX_SIZE);

    return;
}

void write_xtc_frame(t_fileio *traj,unsigned long int step,int natoms,atom *at,rvec *x){


    matrix box;
    int i,j;
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            if(j==i){
                box[i][j] = BOX_SIZE;
            }
            else{
                box[i][j] = 0.0;
            }
        }
    }

    write_xtc(traj,natoms,step,(real) step,box,x,XTC_PREC);

    return;
}


void enforce_constraint(rvec *x,atom *at,molecule mol,rvec ini){

    rvec pos,shift;
    get_bead_center(mol,x,at,pos);
    rvec_sub(ini,pos,shift);
    translate_mol(mol,shift,&x);

    return;
}

int accept_or_reject(double delta,double inv_kT){

    if(delta<=0){
        return 1;
    }
    else{
        double boltz=exp(-delta*inv_kT);
        double r=((double)rand())/((double) RAND_MAX);

        if(boltz>r){
            return 1;
        }
        else{
            return 0;
        }
    }
}

void do_move(params *p,int move,molecule mol,dihedmap dhm,atom *at,int **bond,rvec *old,rvec *new){

    double r=((double)rand())/((double) RAND_MAX)*2.0-1.0;
    int ri,i;
    int *list;
    int nrot=0;
    int active;
    rvec origin,axis;
    clear_rvec(origin);
    clear_rvec(axis);

    
        if(move == 0 || dhm.nactive == 0){
             ri=rand()%3;
             axis[ri]=1.0;

             copy_rvec(old[mol.first],origin);

             list=(int *) calloc(mol.nat,sizeof(int));
             for(i=0;i<mol.nat;i++){
                 list[i]=mol.first+i;
             }

             r*=p->max_rot[move];

             rotate(mol.nat,list,axis,r,origin,old,new);
        }
        else if(move == 1){

             prepare_mapped_dih_rot(dhm,mol,old,&axis,&nrot,&list,&active);
             
             copy_rvec(old[active],origin);

             r*=p->max_rot[move];

             rotate(nrot,list,axis,r,origin,old,new);
         }
    

    free(list);

    return;
}


void do_mc(t_topology *top,params *p,int natoms,atom *at,
        rvec *x,int **bonds,int nmol,molecule *mol,double ***nonb,double ***nonbf,
        dihedmap *dhm,int bAppend,char *fn_trjout,int trj_type){

    int i,j;
    unsigned long int step=0,nstep=p->nstep,acc_count=0;
    double bond_pot[nmol],tot_bond,nonb_pot,total_pot;
    double new_bond[nmol],new_tb,new_nonb,new_tot;
    double cforce=0;
    double pot_avg=0,bd_avg=0,nb_avg=0,cf_avg=0,inv_n;
    double inv_kT = p->beta;
    int ri;
    double r;
    int move; // 0 is rotation, 1 is torsion
    int trajfrq=p->trajfreq;
    int logfrq=p->logfreq;
    int flogfrq=p->flogfreq;

    rvec ini_pos[nmol];
    rvec pos[2];
    rvec distance;


    rvec new[natoms];
    rvec f[natoms];
    char iomode[2];

    FILE *log;
    FILE *gro_traj;
    t_fileio *xtc_traj;
    FILE *flog;
    unsigned long int start=1;
    
    if(bAppend == 0){
        iomode[0] = 'w';
        iomode[1] = '\0';
        for(i=0;i<natoms;i++){
            copy_rvec(x[i],new[i]);
        }
    }
    else if(bAppend == 1){
        iomode[0] = 'a';
        iomode[1] = '\0';
        get_continuation(trj_type,"mclog",fn_trjout,"flog",x,natoms,&step);
        start = step+1;
    }
    else if(bAppend == 3){
        iomode[0] = 'a';
        iomode[1] = '\0';
        get_continuation(trj_type,"mclog",fn_trjout,"flog",NULL,natoms,&step);
        start = step+1;
    }

    do_nbforce(p,f,natoms,nonbf,x);

    init_mc(&(top->idef),p,bond_pot,&tot_bond,&nonb_pot,ini_pos,natoms,at,x,nmol,mol,nonb);

    rvec_sub(ini_pos[0],ini_pos[1],distance);

    cforce=constraint_force(mol[0],distance,f);
    cforce-=constraint_force(mol[1],distance,f);
    cforce*=0.5;

    fprintf(stderr,"Initial:\nNonb_Pot: %lg  Bond[1]: %lg Bond[2]: %lg\nPos 1: %f %f %f\nPos 2: %f %f %f Distance: %f\n",
           nonb_pot,bond_pot[0],bond_pot[1],ini_pos[0][0],ini_pos[0][1],ini_pos[0][2],ini_pos[1][0],ini_pos[1][1],ini_pos[1][2],norm(distance));

    fprintf(stderr,"\nFaites vos jeux!\n\n");

    
    log=fopen("mclog",iomode);
    if(trj_type == 0){
        gro_traj=fopen(fn_trjout,iomode);
    }
    else if(trj_type == 1){
        xtc_traj = open_xtc(fn_trjout,iomode);
    }
    flog=fopen("flog",iomode);
    
    total_pot=tot_bond+nonb_pot;
    pot_avg=total_pot;
    bd_avg=tot_bond;
    nb_avg=nonb_pot;
    cf_avg=cforce;


    if(bAppend == 0){
        if(trj_type == 0){
            fprintf(log,"#step  Pot-Energy[inst/avg] Nonb[inst/avg] Bond[inst/avg] Dist\n");
            fprintf(log,"%li %lg %lg %lg %lg %lg %lg %f\n",step,total_pot,pot_avg,nonb_pot,nb_avg,tot_bond,bd_avg,norm(distance));

            write_gro_frame(gro_traj,step,natoms,at,x);
        }
        else if(trj_type == 1){
            write_xtc_frame(xtc_traj,step,natoms,at,x);
        }
    }
    
    fprintf(flog,"%i %lg %lg\n",0,cforce,cf_avg);


    for(step=start;step<=nstep;step++){

        for(i=0;i<natoms;i++){
            copy_rvec(x[i],new[i]);
        }

        new_tot=0;
        new_tb=0;

        r=((double)rand())/((double) RAND_MAX);
        if(r<p->p_rot){
            move=0;
        }
        else{
            move=1;
        }
        
        ri=rand()%nmol;

        do_move(p,move,mol[ri],dhm[ri],at,bonds,x,new);

        enforce_constraint(new,at,mol[ri],ini_pos[ri]);

        new_nonb = do_nonb(p,natoms,nonb,new);

        new_tot += new_nonb;

        for(i=0;i<nmol;i++){
            new_bond[i]=do_bonded(mol[i].first,mol[i].first+mol[i].nat,new,&(top->idef));
            new_tot += new_bond[i];
            new_tb += new_bond[i];
        }

        if(accept_or_reject(new_tot-total_pot,inv_kT)){

            acc_count+=1;
            total_pot = new_tot;
            nonb_pot = new_nonb;
            for(i=0;i<nmol;i++){
                bond_pot[i]=new_bond[i];
            }
            tot_bond = new_tb;
            for(i=0;i<natoms;i++){
                copy_rvec(new[i],x[i]);
            }
        }

        get_bead_center(mol[0],x,at,pos[0]);
        get_bead_center(mol[1],x,at,pos[1]);
        rvec_sub(pos[0],pos[1],distance);

        inv_n = 1.0 /( (double) (step+1) );
        pot_avg+=total_pot;
        bd_avg+=tot_bond;
        nb_avg+=nonb_pot;

        if(!(step%flogfrq)  || step==nstep){
            do_nbforce(p,f,natoms,nonbf,x);
            cforce=constraint_force(mol[0],distance,f);
            cforce-=constraint_force(mol[1],distance,f);
            cforce*=0.5;
            cf_avg+=cforce;

            fprintf(flog,"%li %lg %lg\n",step,cforce,cf_avg*inv_n);
        }

        if(!(step%logfrq) || step==nstep){
            fprintf(log,"%li %lg %lg %lg %lg %lg %lg %f\n",step,total_pot,pot_avg*inv_n,nonb_pot,nb_avg*inv_n,tot_bond,bd_avg*inv_n,norm(distance));
        }
        if(!(step%trajfrq) || step==nstep){
                if(trj_type == 0){
                    write_gro_frame(gro_traj,step,natoms,at,x);
                }
                else if(trj_type == 1){
                    write_xtc_frame(xtc_traj,step,natoms,at,x);                    
                }

        }
    }

    FILE *last = fopen("last.gro","w");
    write_gro_frame(last,step,natoms,at,x);
    fclose(last);

    fclose(log);
    fclose(flog);

    if(trj_type == 0){
        fclose(gro_traj);
    }
    else if(trj_type == 1){
        close_xtc(xtc_traj);
    }

    fprintf(stderr,"Performed %li steps. Accepted %li (%f %%)\n",nstep,acc_count,100.0*((float) acc_count)/((float) nstep));

    return;
}

