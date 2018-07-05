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
#include <string.h>
#include <math.h>


#include "tpxio.h"
#include "crw_types.h"


void read_param_ir(params *p,t_inputrec *ir){

    p->nstep = ir->userint1;
    p->kT = ir->opts.ref_t[0]*8.31447*0.001;
    p->beta = 1/p->kT;
    p->trajfreq = ir->userint2;
    p->logfreq = ir->userint3;
    p->flogfreq = ir->userint4;
    p->p_rot = ir->userreal1;
    p->r_c[0] = ir->rcoulomb;
    p->r_c[1] = ir->rvdw;
    p->max_rot[0] = ir->userreal2 / 180.0 * M_PI;
    p->max_rot[1] = ir->userreal3 / 180.0 * M_PI;

    return;

}


void parse_cmdline(int ac,char *av[],char **fn_confin,char **fn_trjout,int *trj_type,int *bAppend,char **fn_tprin,int *excl_mode){
// parses cmd-line returns name of start configuration, and append status (0: no;1: yes;2: append but start from new configuration)

    int i=1;
    *bAppend=0;

    while(i<ac){

        if(strcmp(av[i],"-c") == 0){
            i++;
            *fn_confin = (char *) calloc(strlen(av[i])+1,sizeof(char));
            strcpy(*fn_confin,av[i]); //No further check for file format. Simulation will not run without consistent input anyway
            i++;
        }
	    else if(strcmp(av[i],"-o") == 0){
            i++;
            *fn_trjout = (char *) calloc(strlen(av[i])+1,sizeof(char));
            strcpy(*fn_trjout,av[i]);
            if ( strcmp(  &((*fn_trjout)[strlen(*fn_trjout)-3]) , "gro" ) == 0){
                *trj_type = 0;
            }
            else if ( strcmp(  &((*fn_trjout)[strlen(*fn_trjout)-3]) , "xtc" ) == 0  ){
                *trj_type = 1;
            }
            else {
                fprintf(stderr,"WARNING: Unknown trajectory file type. Assuming .gro!\n");
                *trj_type = 0;
            }
            i++;
	    }
        else if(strcmp(av[i],"-append") == 0){
            *bAppend = 1;
            i++;
        }
        else if(strcmp(av[i],"-append_withconf") == 0){
            *bAppend = 2;
            i++;
        }
        else if(strcmp(av[i],"-excl_q") == 0){
            *excl_mode = 1;
            i++;
        }
        else if(strcmp(av[i],"-excl_vdw") == 0){
            *excl_mode = 2;
            i++;
        }
        else if(strcmp(av[i],"-s") == 0){
            i++;
            *fn_tprin = (char *) calloc(strlen(av[i])+1,sizeof(char));
            strcpy(*fn_tprin,av[i]); 
            i++;
        }
        else{
            fprintf(stderr,"Error: Unknown command line option %s\n",av[i]);
            exit(3);
            i++;
        }

    }

    return;
}
