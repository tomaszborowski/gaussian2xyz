#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:09:10 2021

This script reads Gaussian log file for geometry 1-D scan or IRC calculations
and outputs xyz file with optimized geometries for points along the scanned
coordinate or IRC.

@author: Tomasz Borowski
"""
import sys
import matplotlib.pyplot as plt
from extract_geoms_aux import log_read_geo, log_read_step_number_line
from extract_geoms_aux import log_irc_or_scan, log_read_scf, is_geom_converged
from extract_geoms_aux import log_is_ONIOM, log_read_oniom_e


### ---------------------------------------------------------------------- ###
### Seting the file names                                                  ###
inp_file_name = sys.argv[1]
fig_file_name = inp_file_name + ".png"

### ---------------------------------------------------------------------- ###
### test cases
#inp_file_name = './input_examples/oh_h2o.scan.log'


### ---------------------------------------------------------------------- ###
### parsing the log/out file                                               ###
input_f = open(inp_file_name, 'r')

RUN_TYPE = log_irc_or_scan(input_f)
ONIOM = log_is_ONIOM(input_f)

if RUN_TYPE == "SCAN":
    scan_geometries = []
    temp_geo = None
    while temp_geo != "EOF":
        temp_geo = log_read_geo(input_f)
        if temp_geo != "EOF":
            if ONIOM:
                e = log_read_oniom_e(input_f)
            else:
                e = log_read_scf(input_f)
            slt = log_read_step_number_line(input_f)
            if slt:
                last_step = (slt.step_nr == slt.step_max)
                if is_geom_converged(input_f) or last_step:
                    temp_geo.set_scfenergy(e)
                    temp_geo.set_scan_point(slt.scan_point)
                    scan_geometries.append(temp_geo)
    
elif RUN_TYPE == "IRC":
    pass


input_f.close()


### ---------------------------------------------------------------------- ###
### generating the output                                                  ###

if RUN_TYPE == "SCAN":
    points = []
    energie = []
    for geo in scan_geometries:
        geo.print_xyz()
        points.append( geo.get_scan_point() )
        energie.append( geo.get_scfenergy() )
        
    plt.figure(1)
    plt.plot(points, energie, 'go--', linewidth=1, markersize=6)
    plt.grid(b=True, which='major', axis='both')
    plt.xticks(ticks=points)
    plt.xlabel('Scan point number')
    plt.ylabel('SCF/ONIOM E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)
    

elif RUN_TYPE == "IRC":
    pass