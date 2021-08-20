#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:09:10 2021

This script reads Gaussian log file and outputs xyz with optimized geometry/ies.

For geometry 1-D scan or IRC calculations optimized (or last) geometries 
for points along the scanned coordinate or IRC are output. A png file with a plot
of ONIOM or SCF energies along the profile is generated.

The script expects 2 or 3 arguments: 
    #1 log-file-name, 
    #2 type of extraction - one from amoung: scan, irc, all, last, 
    #3 (for IRC) file name with SP/FREQ calculations for the TS from which IRC calculations started

@author: Tomasz Borowski
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from extract_geoms_aux import log_read_geo, log_read_step_number_line
from extract_geoms_aux import log_irc_or_scan, log_read_scf, is_geom_converged
from extract_geoms_aux import log_is_ONIOM, log_read_oniom_e, is_irc_converged
from extract_geoms_aux import log_read_irc_data

LEGIT_RUN_TYPE = ["SCAN", "IRC", "ALL", "LAST"]

def read_geo_scf_oniom_e(input_f):
    """ from Gaussian output file read geometry and its SCF and ONIOM energies 
    INPUT: input_f - file object (Gaussian output file)
    OUTPUT: temp_geo - Geometry object or "EOF" string  """
    temp_geo = log_read_geo(input_f)
    scf_e = None
    oniom_e = None
    if temp_geo != "EOF":
        scf_e = log_read_scf(input_f)
        temp_geo.set_scf_energy(scf_e)
        if ONIOM:
            oniom_e = log_read_oniom_e(input_f)
            temp_geo.set_oniom_energy(oniom_e)
    return temp_geo


RUN_TYPE = None

### ---------------------------------------------------------------------- ###
### Seting the file names                                                  ###
inp_file_name = sys.argv[1]
fig_file_name = inp_file_name + ".png"

if sys.argv[2]:
    flag_read = sys.argv[2]
    if isinstance(flag_read, str):
        flag_read = flag_read.upper()
        if flag_read in LEGIT_RUN_TYPE:
            RUN_TYPE = flag_read

# optionally to read TS geometry and energy (IRC point 0)
irc_ts_file_name = None            
if sys.argv[3]:
    irc_ts_file_name = sys.argv[3]


### ---------------------------------------------------------------------- ###
### test cases
#inp_file_name = './input_examples/oh_h2o.scan.log'
#inp_file_name = './input_examples/2x_scan.log'
#inp_file_name = './input_examples/h2o_opt.log'
#inp_file_name = './input_examples/h2o_sp.log'
# inp_file_name = './input_examples/oh_h2o.irc.log'

# RUN_TYPE = "IRC"


### ---------------------------------------------------------------------- ###
### Fig file names                                                         ###
fig_file_name = inp_file_name + ".png"



### ---------------------------------------------------------------------- ###
### parsing the log/out file                                               ###
input_f = open(inp_file_name, 'r')

if RUN_TYPE not in LEGIT_RUN_TYPE:
    RUN_TYPE = log_irc_or_scan(input_f)
    
ONIOM = log_is_ONIOM(input_f)


if RUN_TYPE == "SCAN":
    scan_geometries = []
    temp_geo = None
    while temp_geo != "EOF":
        temp_geo = read_geo_scf_oniom_e(input_f)                
        slt = log_read_step_number_line(input_f)
        if slt:
            temp_geo.set_scan_point( slt.scan_point )
            geom_conv = is_geom_converged(input_f)
            temp_geo.set_geom_converged( geom_conv )
            last_step = (slt.step_nr == slt.step_max)            
            if geom_conv or last_step:
                temp_geo.set_in_scan(True)
                scan_geometries.append(temp_geo)

    
elif RUN_TYPE == "IRC":
    irc_geometries = [] 
    temp_geo = None
    while temp_geo != "EOF":
        temp_geo = read_geo_scf_oniom_e(input_f)
        irc_conv = is_irc_converged(input_f)
        if irc_conv:
            irc_pt = log_read_irc_data(input_f)
            if irc_pt:
                temp_geo.set_irc_path_number( irc_pt.path_nr )
                temp_geo.set_irc_point_number( irc_pt.point_nr )
                temp_geo.set_irc_net_reaction_coordinate( irc_pt.net_reaction_coord )
                temp_geo.set_in_irc(True)
                irc_geometries.append(temp_geo)
            

if RUN_TYPE == "ALL":
    seq_nr = []
    energie = []
    i = 1
    temp_geo = None
    while temp_geo != "EOF":
        temp_geo = read_geo_scf_oniom_e(input_f)
        if temp_geo != "EOF" and temp_geo.get_scf_energy():                
            temp_geo.print_xyz()
            seq_nr.append( i )
            i += 1
            energie.append( temp_geo.get_scf_energy() )


if RUN_TYPE =="LAST":
    temp_geo = None
    while temp_geo != "EOF":
        if temp_geo == None:
            prev_temp_geo = None
        elif temp_geo.get_scf_energy():
            prev_temp_geo = temp_geo
        temp_geo = read_geo_scf_oniom_e(input_f)
        if temp_geo == "EOF":
            if prev_temp_geo != "EOF" and prev_temp_geo != None:
                prev_temp_geo.print_xyz()


input_f.close()


### ---------------------------------------------------------------------- ###
### optionally for IRC read the TS from a separate file                    ###

if RUN_TYPE == "IRC" and irc_ts_file_name:    
    irc_ts_f = open(irc_ts_file_name, 'r')
    temp_geo = read_geo_scf_oniom_e(irc_ts_f)
    temp_geo.set_irc_path_number( 1 )
    temp_geo.set_irc_point_number( 0 )
    temp_geo.set_irc_net_reaction_coordinate( 0.0 )
    temp_geo.set_in_irc(True)
    irc_geometries.append(temp_geo)    

irc_ts_f.close()

### ---------------------------------------------------------------------- ###
### generating the output                                                  ###

if RUN_TYPE == "SCAN":
    seq_nr = []
    points = []
    energie = []
    i = 1
    for geo in scan_geometries:
        geo.print_xyz()
        seq_nr.append( i )
        i += 1
        points.append( geo.get_scan_point() )
        if ONIOM:
            energie.append( geo.get_oniom_energy() )
        else:
            energie.append( geo.get_scf_energy() )
        
        
    plt.figure(1)
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True)) 
    plt.plot(seq_nr, energie, 'go--', linewidth=1, markersize=6)   
    plt.grid(b=True, which='major', axis='both')
    plt.xlabel('Scan point number')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)
    

elif RUN_TYPE == "IRC":
    irc_net_coord = []
    energie = []
    for geo in irc_geometries:
        if geo.get_irc_path_number() == 2:
            net_coord = -1.0 * geo.get_irc_net_reaction_coordinate() # assumes path nr 2 is for negative rection coordinate
            geo.set_irc_net_reaction_coordinate( net_coord )
    
    irc_geometries.sort( key=lambda x: x.get_irc_net_reaction_coordinate() )
    
    for geo in irc_geometries:
        geo.print_xyz()
        irc_net_coord.append( geo.get_irc_net_reaction_coordinate() )
        if ONIOM:
            energie.append( geo.get_oniom_energy() )
        else:
            energie.append( geo.get_scf_energy() )
 
    plt.figure(1)
    plt.plot(irc_net_coord, energie, 'go--', linewidth=1, markersize=6)   
    plt.grid(b=True, which='major', axis='both')
    plt.xlabel('IRC net coordinate')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)
       

elif RUN_TYPE == "ALL":
    plt.figure(1)
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.plot(seq_nr, energie, 'go--', linewidth=1, markersize=6)
    plt.grid(b=True, which='major', axis='both')
    plt.xlabel('Geometry number')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)   

