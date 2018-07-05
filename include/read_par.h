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

#ifndef READ_PAR_H
#define READ_PAR_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int map_str(char *s);

void read_param_ir(params *p,t_inputrec *ir);
void read_param(params *p);
void parse_cmdline(int ac,char *av[],char **fn_confin,char **fn_trjout,int *trj_type,int *bAppend,char **fn_tprin,int *excl_mode);

#endif
