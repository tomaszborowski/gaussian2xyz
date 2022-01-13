#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 16:09:10 2021

This script reads Gaussian log file and outputs xyz with optimized geometry/ies.

For geometry 1-D scan or IRC calculations optimized (or last) geometries 
for points along the scanned coordinate or IRC are output. A png file with a plot
of ONIOM or SCF energies along the profile is generated.

The script expects 2 or 3 arguments: 
    #1 log-file-name or, for irc_f, name of a file specifying filenames and directions of irc calc.
    #2 type of extraction - one from amoung: scan, irc, irc_f, all, last, 
    #3 (for IRC) file name with SP/FREQ calculations for the TS from which IRC calculations started

@author: Tomasz Borowski, Zuzanna Wojdyła

branch: zuza
last modification: 13.01.2022
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from extract_geoms_aux import log_read_geo, log_read_step_number_line
from extract_geoms_aux import log_irc_or_scan, log_read_scf, is_geom_converged
from extract_geoms_aux import log_is_ONIOM, log_read_oniom_e, is_irc_converged
from extract_geoms_aux import log_read_irc_data

LEGIT_RUN_TYPE = ["SCAN", "IRC", "IRC_F", "ALL", "LAST"]
LEGIT_DIRECTIONS = ["REVERSE", "TS", "FORWARD"]

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
energy_file_name = inp_file_name + ".dat"

if sys.argv[2]:
    flag_read = sys.argv[2]
    if isinstance(flag_read, str):
        flag_read = flag_read.upper()
        if flag_read in LEGIT_RUN_TYPE:
            RUN_TYPE = flag_read

# optionally to read TS geometry and energy (IRC point 0)
irc_ts_file_name = None            
if len(sys.argv) > 3:
    irc_ts_file_name = sys.argv[3]

#print(sys.argv)

### ---------------------------------------------------------------------- ###
### test cases
#inp_file_name = './input_examples/oh_h2o.scan.log'
#inp_file_name = './input_examples/2x_scan.log'
#inp_file_name = './input_examples/h2o_opt.log'
#inp_file_name = './input_examples/h2o_sp.log'
#inp_file_name = './input_examples/oh_h2o.irc.log'

# irc_ts_file_name = None
# inp_file_name = './input_examples/oh_h2o.irc.info'
# RUN_TYPE = "IRC_F"

# irc_ts_file_name = None
# inp_file_name = './input_examples/oh_h2o.irc_1b.log'
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
    irc_last_point = None
    temp_geo = None
    while temp_geo != "EOF":
        temp_geo = read_geo_scf_oniom_e(input_f)
        irc_conv = is_irc_converged(input_f)
        if irc_conv:
            irc_pt = log_read_irc_data(input_f, irc_last_point)
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
### case IRC_F                                                             ###
irc_file_names = []
irc_directions = [] 
irc_last = [] 
                                                 
if RUN_TYPE == "IRC_F":   
    while True:
        a = input_f.readline()
        if not a:
            break
        a_split = a.split()
        irc_file_names.append(a_split[0])
        direction_read = a_split[1].upper()
        if len(a.split())==3:
            irc_last.append(a_split[2])
        else:
            irc_last.append(None)
        if direction_read in LEGIT_DIRECTIONS:
            irc_directions.append(direction_read)
input_f.close()

if RUN_TYPE == "IRC_F":
    nr_irc_files = len(irc_file_names)
    nr_irc_directions = len(irc_directions)
    
    if nr_irc_files != nr_irc_directions:
        print("\n Huston, we've got a problem with a file specifying irc files\n")
        sys.exit("I am exiting")
    else:
        irc_geometries = [] 
        react_coord_offset = 0.0
        last_reaction_coord = 0.0
        point_offset = 0
        last_point_nr = 0
        prev_path = None
        for f_name, irc_dir, irc_last_point in zip(irc_file_names, irc_directions, irc_last):
            temp_geo = None
            if irc_dir == 'REVERSE':
                path_nr = 2
            elif irc_dir == 'FORWARD':
                path_nr = 1
            if path_nr != prev_path:
                prev_path = path_nr
                react_coord_offset = 0.0
                point_offset = 0
            else:
                react_coord_offset = last_reaction_coord
                point_offset = last_point_nr
            input_f = open(f_name, 'r')
            if irc_dir == 'TS':
                temp_geo = read_geo_scf_oniom_e(input_f)
                temp_geo.set_irc_path_number( 1 )
                temp_geo.set_irc_point_number( 0 )
                temp_geo.set_irc_net_reaction_coordinate( 0.0 )
                temp_geo.set_in_irc(True)
                irc_geometries.append(temp_geo)
            else:    
                while temp_geo != "EOF":
                    temp_geo = read_geo_scf_oniom_e(input_f)
                    irc_conv = is_irc_converged(input_f)
                    if irc_conv:
                        irc_pt = log_read_irc_data(input_f,irc_last_point)
                        if irc_pt:
                            temp_geo.set_irc_path_number( path_nr )                       
                            temp_geo.set_irc_point_number( irc_pt.point_nr + point_offset )
                            temp_geo.set_irc_net_reaction_coordinate( irc_pt.net_reaction_coord  + react_coord_offset)
                            temp_geo.set_in_irc(True)
                            irc_geometries.append(temp_geo)
                            last_reaction_coord = irc_pt.net_reaction_coord + react_coord_offset
                            last_point_nr = irc_pt.point_nr + point_offset
            input_f.close()
        
        
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
#    plt.grid(b=True, which='major', axis='both')
    plt.grid(visible=True, which='major', axis='both')
    plt.xlabel('Scan point number')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)
    

elif RUN_TYPE == "IRC" or RUN_TYPE == "IRC_F":
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
    #normalise irc coord
    min_val=min(irc_net_coord)
    max_val=max(irc_net_coord)
    irc_normalised_coord = [(x-min_val)/(max_val-min_val) for x in irc_net_coord]

    plt.figure(1)
    ###plot for normalised IRC coord or one in amu**0.5 Bohr (?) ###
    #plt.plot(irc_normalised_coord, energie, 'go--', linewidth=1, markersize=6)   
    plt.plot(irc_net_coord, energie, 'go--', linewidth=1, markersize=6)   
#    plt.grid(b=True, which='major', axis='both')
    plt.grid(visible=True, which='major', axis='both')
    plt.xlabel('IRC net coordinate')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)
       
    ### write energy and coords to file ###
    energy_file=open(energy_file_name, 'w')
    for i in range(len(irc_net_coord)):
        line=str("{:.6f}".format(irc_net_coord[i]))+"\t"+str("{:.6f}".format(irc_normalised_coord[i]))+"\t"+str(energie[i])+"\n"
        energy_file.write(line)
    energy_file.close()

elif RUN_TYPE == "ALL":
    plt.figure(1)
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.plot(seq_nr, energie, 'go--', linewidth=1, markersize=6)
#    plt.grid(b=True, which='major', axis='both')
    plt.grid(visible=True, which='major', axis='both')
    plt.xlabel('Geometry number')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)   

