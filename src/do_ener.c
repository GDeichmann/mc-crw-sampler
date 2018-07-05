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


#include "do_ener.h"

double do_bonded(int min,int max,rvec *x,t_idef *id){ 
//calculates all (implemented) bonded pot. energies in idef for one molecule: 
//for min and max enter the first atom of a molecule and the first atom of the next molecule.
//x is the _GLOBAL_ coordinate array                                                      
    
    double pot=0.0;
    int at_nr,par_nr,type,i,j,k,l;
    int block; //blocks for bond,angle and dihedrals
    int at[4]; //there will be max. 4 atoms and they wil be named a b c d with distance vectors ab bc and cd
    rvec ab,bc,cd;
    double d_ab,d_bc,d_cd; //distances
    double a_abc,a_bcd; // angles
    rvec p_abc,p_bcd; //plane normal vectors
    double dih,cos_dih; //dihedral
    double rb_par[6]; //for fourier to rb conversion
    t_iparams ip;
    const double deg_to_rad = M_PI/180.0;

    double pot_prev=0.0;

    for(i=0;i<F_NRE;i++){//For complete list see GMXINCLUDE/gromacs/types/idef.h

        if(i>=F_BONDS  && i<=F_RESTRBONDS){
            at_nr=2;
            block=1;
        }
        else if(i>=F_ANGLES  && i<=F_TABANGLES){
            at_nr=3;
            block=2;
        }
        else if(i>=F_PDIHS && i<=F_TABDIHS){
            at_nr=4;
            block=3;
        }
        else{ // The rest are nb params and exotic functional types do not use or implement at will
            at_nr=0;
            block=0;
        }

        if(at_nr && id->il[i].nr){

            for(j=0;j<id->il[i].nr;){

                type=id->il[i].iatoms[j];
                at[0]=id->il[i].iatoms[j+1];
                if(at[0] >= min && at[0] < max){ //Assumes no bonded Dof between molecules; 
                                                                     //Could get rid of this but then we do not distiguish between molecules
                    ip = id->iparams[type];
                    for(k=1;k<at_nr;k++){
                        at[k]=id->il[i].iatoms[j+k+1];
                    }
                    
                    if(block == 1){
                        rvec_sub(x[at[0]],x[at[1]],ab);
                        d_ab = norm(ab);

                        switch(i){
                            case F_BONDS:
                            case F_HARMONIC:
                                pot += 0.5*ip.harmonic.krA*(d_ab-ip.harmonic.rA)*(d_ab-ip.harmonic.rA);
                                break;
                            case F_G96BONDS:
                                pot += 0.25*ip.harmonic.krA*(d_ab*d_ab-ip.harmonic.rA)*(d_ab*d_ab-ip.harmonic.rA) ;
                                break;
                            case F_MORSE:
                                pot += ip.morse.cbA*(1-exp(-ip.morse.betaA*(d_ab-ip.morse.b0A)))*(1-exp(-ip.morse.betaA*(d_ab-ip.morse.b0A)));
                                break;
                            case F_CUBICBONDS:
                                pot += ip.cubic.kb*(d_ab-ip.cubic.b0)*(d_ab-ip.cubic.b0) + 
                                    ip.cubic.kb*ip.cubic.kcub*(d_ab-ip.cubic.b0)*(d_ab-ip.cubic.b0)*(d_ab-ip.cubic.b0);
                                break;
                            case F_FENEBONDS:
                            case F_TABBONDS:
                            case F_TABBONDSNC:
                            case F_CONNBONDS:                            
                            default:
                                fprintf(stderr,"Exotic bond type no. %i not implemented\nFor more look up hash table in GMXINCLUDE/gromacs/types/idef.h\n",i);
                        }                                

                    }
                    else if(block == 2){
                        rvec_sub(x[at[0]],x[at[1]],ab);
                        rvec_sub(x[at[1]],x[at[2]],bc);
                        
                        switch(i){
                            case F_ANGLES:
                                a_abc = M_PI - gmx_angle(ab,bc);
                                pot += 0.5*ip.harmonic.krA*( a_abc - (ip.harmonic.rA)*deg_to_rad  )*( a_abc - (ip.harmonic.rA)*deg_to_rad );
                                break;
                            case F_G96ANGLES:
                                a_abc = cos_angle(ab,bc);
                                pot += 0.5*ip.harmonic.krA*( (a_abc + ip.harmonic.rA) * (a_abc + ip.harmonic.rA)  );
                                break;
                            case F_LINEAR_ANGLES:
                            case F_CROSS_BOND_BONDS:
                            case F_CROSS_BOND_ANGLES:
                            case F_UREY_BRADLEY:
                            case F_QUARTIC_ANGLES:
                            case F_TABANGLES:
                            default:
                                fprintf(stderr,"Exotic bond type no. %i not implemented\nFor more look up hash table in GMXINCLUDE/gromacs/types/idef.h\n",i);
                        }       


                    }
                    else if(block == 3){
                        rvec_sub(x[at[0]],x[at[1]],ab);
                        rvec_sub(x[at[1]],x[at[2]],bc);
                        rvec_sub(x[at[2]],x[at[3]],cd);
                        
                        cprod(ab,bc,p_abc);
                        cprod(bc,cd,p_bcd);
                        
                        switch(i){
                            case F_PDIHS:
                                dih = gmx_angle(p_abc,p_bcd);
                                pot += ip.pdihs.cpA*(1.0+cos((float)(ip.pdihs.mult)*dih - ip.pdihs.phiA*deg_to_rad));
                                break;
                            case F_RBDIHS:
                            case F_FOURDIHS:
                                cos_dih=-1.0*cos_angle(p_abc,p_bcd);//conversion between conventions takes place here; takes advantage of symmetry of cosine function
                                for(l=0;l<6;l++){
                                    double inc=ip.rbdihs.rbcA[l]*pow(cos_dih,l);
                                    pot+=inc;
                                    }
                                break;
                            case F_IDIHS:
                                dih = gmx_angle(p_abc,p_bcd);
                                pot += 0.5*ip.harmonic.krA*( dih - (ip.harmonic.rA)*deg_to_rad  )*( dih - (ip.harmonic.rA)*deg_to_rad );
                                break;
                            case F_PIDIHS:
                            case F_TABDIHS:
                            default:
                                fprintf(stderr,"Exotic bond type no. %i not implemented\nFor more look up hash table in GMXINCLUDE/gromacs/types/idef.h\n",i);
                        }
                    }

                }

                j+=at_nr+1;

            }

        }

    }
    


    return pot;
}


double do_nonb(params *p,int natoms,double ***nb,rvec *x){

    double coul=0,disp=0,rep=0;

    int i,j;
    rvec r;
    real dist_sq,inv_dist;
    double inv_6;

    for(i=0;i<natoms;i++){
        for(j=i+1;j<natoms;j++){
        //for(j=0;j<i;j++){
            if(i!=j){
                rvec_sub(x[i],x[j],r);
                dist_sq = norm2(r);
                inv_dist = gmx_invsqrt(dist_sq);
                inv_6 = 1.0 / (dist_sq*dist_sq*dist_sq);

                if(p->r_c[0] < 0 || dist_sq <= (p->r_c[0]*p->r_c[0])){
                    coul+=nb[0][i][j]*inv_dist;
                }
                if(p->r_c[1] < 0 || dist_sq <= (p->r_c[1]*p->r_c[1])){
                    disp+=nb[1][i][j]*inv_6;
                    rep+=nb[2][i][j]*inv_6*inv_6;
                }
                
            }
        }
    }
    return coul-disp+rep;
}

