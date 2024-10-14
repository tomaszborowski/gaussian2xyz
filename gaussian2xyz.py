#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script reads Gaussian log file and outputs xyz geometry/ies.

For geometry 1-D scan or IRC calculations optimized (or last) geometries 
for points along the scanned coordinate or IRC are output. A png file with a plot
of ONIOM or SCF energies along the profile is generated.

The script expects 2 or 3 arguments: 
    #1 log-file-name or, for irc_f, name of a file specifying filenames and directions of irc calc.
    #2 type of extraction - one from amoung: scan, irc, irc_f, all, last, nr, lasts, nrs
    #3 (for IRC) file name with SP/FREQ calculations for the TS from which IRC calculations started
       or for NR 1-based number of the structure to be extracted

@authors: Tomasz Borowski, Zuzanna Wojdyła
modification: 17.10.2022, 18.05.2023, 19.05.2023, 23.05.2023, 31.10.2023, 2.11.2023 
last modification: 14.10.2024
"""

import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from extract_geoms_aux import log_read_step_number_line
from extract_geoms_aux import is_geom_converged, scan_file
from extract_geoms_aux import log_is_ONIOM, is_irc_converged
from extract_geoms_aux import log_read_irc_data, log_is_MM, print_help
from extract_geoms_aux import log_read_inp_geo, log_read_mulliken
from extract_geoms_aux import read_geo_scf_oniom_e
from extract_geoms_aux import read_geo_mm_e, find_first_infile

LEGIT_RUN_TYPE = ["SCAN", "IRC", "IRC_F", "ALL", "LAST", "NR", "LASTS", "NRS"]
LEGIT_DIRECTIONS = ["REVERSE", "TS", "FORWARD"]

           
### ---------------------------------------------------------------------- ###
RUN_TYPE = None
ONIOM = False       # if the log file is from ONIOM calculations
MM = False          # if the log file is from MM calculations
n_atoms = 0         # number of atoms in the system
is_input_orient = False  # if geometry is in Input orientation 
is_standard_orient = False # if geometry is in Standard orientation

### ---------------------------------------------------------------------- ###
### Seting the file names                                                  ###    
sys_argv_len = len(sys.argv)
if sys_argv_len > 1:
    inp_file_name = sys.argv[1]
    fig_file_name = inp_file_name + ".png"
    energy_file_name = inp_file_name + ".dat"
else:
    inp_file_name = None

if inp_file_name == "-h":
    print_help()
    sys.exit(1)

if inp_file_name == None:
    print("log-file-name or, for irc_f, name of a file specifying filenames and directions of irc calc not found \n")
    sys.exit(1)

if sys.argv[2]:
    flag_read = sys.argv[2]
    if isinstance(flag_read, str):
        flag_read = flag_read.upper()
        if flag_read in LEGIT_RUN_TYPE:
            RUN_TYPE = flag_read
        else:
            print("Flag for extraction type that was provided: " + flag_read + " has not been recognized \n")
            sys.exit(1)


# optionally to read TS geometry and energy (IRC point 0)
# or number of structure ("NR")
irc_ts_file_name = None    
str_number = None        
if len(sys.argv) > 3:
    if RUN_TYPE == "IRC_F":
        irc_ts_file_name = sys.argv[3]
    elif RUN_TYPE == "NR":
        str_number = eval( sys.argv[3] )


### ---------------------------------------------------------------------- ###
### initial parsing the log/out file                                       ###
input_f = open(inp_file_name, 'r')

inp_orient_positions = []
std_orient_positions = []

if RUN_TYPE != "IRC_F":
    if RUN_TYPE =="LAST" or RUN_TYPE =="LASTS" or RUN_TYPE == "NR" or RUN_TYPE == "NRS":
        inp_orient_positions = scan_file(input_f, "Input orientation:")
        if len(inp_orient_positions) > 0:
            is_input_orient = True
        if not is_input_orient:
            std_orient_positions = scan_file(input_f, "Standard orientation:")
            if len(std_orient_positions) > 0:
                is_standard_orient = True
    else:
        inp_orient_positions = find_first_infile(input_f, "Input orientation:")
        if len(inp_orient_positions) > 0:
            is_input_orient = True
        if not is_input_orient:
            std_orient_positions = find_first_infile(input_f, "Standard orientation:")
            if len(std_orient_positions) > 0:
                is_standard_orient = True
                
elif RUN_TYPE == "IRC_F":
    irc_file_names = []
    irc_directions = [] 
    irc_last = [] 
                                                       
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

    file_name_1 = irc_file_names[0] 
    input_f_1 = open(file_name_1, 'r') # using first irc file as a reference
    inp_orient_positions = find_first_infile(input_f_1, "Input orientation:")
    if len(inp_orient_positions) > 0:
        is_input_orient = True
    if not is_input_orient:
        std_orient_positions = find_first_infile(input_f_1, "Standard orientation:")
        if len(std_orient_positions) > 0:
            is_standard_orient = True    

if is_input_orient:
    flag_line = "Input orientation:"
    geometry_positions = inp_orient_positions
elif is_standard_orient:
    flag_line = "Standard orientation:"
    geometry_positions = std_orient_positions
else:
    print("Neither Input nor Standard orientation found!")
    exit(1)

if RUN_TYPE == "LASTS" or RUN_TYPE == "NRS":
    mull_qs_positions = scan_file(input_f, "Mulliken charges and spin densities:")

if RUN_TYPE != "IRC_F":
    ONIOM = log_is_ONIOM(input_f)
    if not ONIOM:
        MM = log_is_MM(input_f)
elif RUN_TYPE == "IRC_F":
    ONIOM = log_is_ONIOM(input_f_1)
    if not ONIOM:
        MM = log_is_MM(input_f_1)    

if MM:
    inp_geo_flag = 'Symbolic Z-Matrix:'
else:
    inp_geo_flag = 'Symbolic Z-matrix:'

if RUN_TYPE != "IRC_F":    
    input_geo = log_read_inp_geo(input_f, inp_geo_flag) # read input geometry (for ONIOM with info about layers and LAH)   
elif RUN_TYPE == "IRC_F":
    input_geo = log_read_inp_geo(input_f_1, inp_geo_flag)
    input_f_1.close()

n_atoms = input_geo.get_n_atoms()
    
### ---------------------------------------------------------------------- ###
### case specific parsing the log/out file                                       ###

if RUN_TYPE == "SCAN":
    scan_geometries = []
    temp_geo = None
    while temp_geo != "EOF": 
        temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)               
        slt = log_read_step_number_line(input_f)
        if slt:
            temp_geo.set_scan_point( slt.scan_point )
            geom_conv = is_geom_converged(input_f)
            temp_geo.set_geom_converged( geom_conv )
            last_step = (slt.step_nr == slt.step_max)            
            if geom_conv or last_step:
                temp_geo.set_in_scan(True)
                scan_geometries.append(temp_geo)
    input_f.close()

    
elif RUN_TYPE == "ALL":
    seq_nr = []
    energie = []
    i = 1
    temp_geo = None
    while temp_geo != "EOF":
        if MM:
            temp_geo = read_geo_mm_e(input_f, flag_line, n_atoms)
        else:
            temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)
        if temp_geo != "EOF" and temp_geo.get_scf_energy():                
            temp_geo.print_xyz()
            seq_nr.append( i )
            i += 1
            if ONIOM:
                energie.append( temp_geo.get_oniom_energy() )
            else:
                energie.append( temp_geo.get_scf_energy() )
    input_f.close()


elif RUN_TYPE =="LAST" or RUN_TYPE =="LASTS" or RUN_TYPE == "NR" or RUN_TYPE == "NRS":
    if MM:
        energy_positions = scan_file(input_f, "Energy per function class:")
    elif ONIOM:
        energy_positions = scan_file(input_f, "ONIOM: extrapolated energy =")
    else:
        energy_positions = scan_file(input_f, "SCF Done:")
    
    if RUN_TYPE =="LAST" or RUN_TYPE =="LASTS":
        wanted_e_pos = energy_positions[-1]
    elif RUN_TYPE == "NR" or RUN_TYPE == "NRS":
        wanted_e_pos = energy_positions[str_number - 1]
    
    if MM: # in MM log energy is reported before the geometry
        jump_pos = wanted_e_pos - 200 # small offset (200) to be sure jump_pos is at lest one line earlier
    else: # in ONIOM / SCF energy is reported after the geometry 
        geometry_positions.reverse()
        for g_pos in geometry_positions:
            if g_pos < wanted_e_pos:
                jump_pos = g_pos - 200
                break
    
    input_f.seek(jump_pos)
    if MM:
        temp_geo = read_geo_mm_e(input_f, flag_line, n_atoms)
    else:
        temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)
    if RUN_TYPE =="LAST" or RUN_TYPE == "NR":
        temp_geo.print_xyz()
    elif RUN_TYPE =="LASTS" or RUN_TYPE == "NRS":
        if ONIOM: 
            mull_qs_positions.reverse() # in ONIOM Mulliken charges and spin pops are reported before ONIOM Energy
            for m_pos in mull_qs_positions:
                if m_pos < wanted_e_pos:
                    jump_pos = m_pos
                    break
            input_f.seek(jump_pos)
            mull_q, mull_s = log_read_mulliken(input_f)
        else:
            for m_pos in mull_qs_positions:
                if m_pos > wanted_e_pos: # in QM (not ONIOM) job logs Mulliken charges and spin pops are reported after SCF Energy
                    jump_pos = m_pos
                    break
            input_f.seek(jump_pos)
            mull_q, mull_s = log_read_mulliken(input_f)            
        if len(mull_q) == 0:
            print("Mulliken charges (and spin populations) not found \n")
            exit(1)
        else:
            temp_geo.set_mulliken_spin_pops(mull_s)
            if ONIOM:
                for at_temp_geo, at_input_geo in zip(temp_geo.get_atoms(), input_geo.get_atoms()):
                    at_temp_geo.set_oniom_layer(at_input_geo.get_oniom_layer())
                    at_temp_geo.set_link_atom_host(at_input_geo.is_link_atom_host())                    
                temp_geo.set_h_lah_atoms()
                qm_plus_lah = temp_geo.get_h_lah_atoms()
                qm_plus_lah.set_oniom_energy(temp_geo.get_oniom_energy())
                qm_plus_lah.print_xyzs()
            else:
                temp_geo.print_xyzs()
    input_f.close()


elif RUN_TYPE == "IRC":
    irc_geometries = [] 
    irc_last_point = None
    temp_geo = None
    while temp_geo != "EOF":
        temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)
        irc_conv = is_irc_converged(input_f)
        if irc_conv:
            irc_pt = log_read_irc_data(input_f, irc_last_point)
            if irc_pt:
                temp_geo.set_irc_path_number( irc_pt.path_nr )
                temp_geo.set_irc_point_number( irc_pt.point_nr )
                temp_geo.set_irc_net_reaction_coordinate( irc_pt.net_reaction_coord )
                temp_geo.set_in_irc(True)
                irc_geometries.append(temp_geo)
    input_f.close()

    if irc_ts_file_name:    # optionally for IRC read the TS from a separate file
        irc_ts_f = open(irc_ts_file_name, 'r')
        temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)
        temp_geo.set_irc_path_number( 1 )
        temp_geo.set_irc_point_number( 0 )
        temp_geo.set_irc_net_reaction_coordinate( 0.0 )
        temp_geo.set_in_irc(True)
        irc_geometries.append(temp_geo)    
    
        irc_ts_f.close()
                                                            

elif RUN_TYPE == "IRC_F":
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
                temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)
                temp_geo.set_irc_path_number( 1 )
                temp_geo.set_irc_point_number( 0 )
                temp_geo.set_irc_net_reaction_coordinate( 0.0 )
                temp_geo.set_in_irc(True)
                irc_geometries.append(temp_geo)
            else:    
                while temp_geo != "EOF":
                    temp_geo = read_geo_scf_oniom_e(input_f, flag_line, n_atoms, ONIOM)
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
    ###plot for normalised IRC coord or one in units used by Gaussian ###
    #plt.plot(irc_normalised_coord, energie, 'go--', linewidth=1, markersize=6)   
    plt.plot(irc_net_coord, energie, 'go--', linewidth=1, markersize=6)   
    plt.grid(visible=True, which='major', axis='both')
    plt.xlabel('IRC net coordinate')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)
       
    ### write energy and coords to a file ###
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
    plt.grid(visible=True, which='major', axis='both')
    plt.xlabel('Geometry number')
    if ONIOM:
        plt.ylabel('ONIOM E [a.u.]')
    elif MM:
        plt.ylabel('MM Energy [a.u.]')
    else:
        plt.ylabel('SCF E [a.u.]')
    plt.savefig(fig_file_name, dpi=300)   

