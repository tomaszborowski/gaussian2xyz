#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
a collection of classes and functions for gaussian2xyz

Authors: Tomasz Borowski, Zuzanna Wojdyla
modification: 7.01.2023, 18.05.2023, 19.05.2023, 23.05.2023, 31.10.2023
last modification: 2.11.2023
"""

import numpy as np
import re
from collections import namedtuple
import fortranformat as ff

Bohr_to_A = 0.52917721
Hartree_to_kcal = 627.509608


# at_num_symbol - dic mapping atomic number to element symbol (up to 86 - Rn)
at_num_symbol = \
    {1: 'H', 2: 'He',
     3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
     11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
     19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
     30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
     37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
     48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
     55: 'Cs', 56: 'Ba', 57: 'La', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
     80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
     58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er',
     69: 'Tm', 70: 'Yb', 71: 'Lu'}

# symbol_at_num - dic mapping symbol to atomic number (up to 86 - Rn)
symbol_at_num = \
    {'H': 1, 'He':2,
     'Li': 3, 'Be':4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
     'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
     'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
     'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
     'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47,
     'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
     'Cs': 55, 'Ba': 56, 'La': 57, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
     'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
     'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68,
     'Tm': 69, 'Yb': 70, 'Lu': 71}

class Atom(object):
    atomic_number = 0
    symbol = ""
    number = 0 # 1-based index of an atom
    coordinates = np.empty((1, 3), dtype=np.float64)  # coordinates in Angstroms
    mulliken_charge = 0.0
    mulliken_spin_pop = 0.0
    oniom_layer = None # None or one from {H, M, L}
    link_atom_host = False


    # the class "constructor"
    def __init__(self, atomic_number, coordinates):
        self.atomic_number = atomic_number
        self.symbol = at_num_symbol[atomic_number]
        self.coordinates = np.array(coordinates, dtype=np.float64)

    def get_symbol(self):
        return self.symbol
    def get_number(self):
        return self.number
    def get_coords(self):
        return self.coordinates
    def get_mulliken_charge(self):
        return self.mulliken_charge
    def get_mulliken_spin_pop(self):
        return self.mulliken_spin_pop
    def get_oniom_layer(self):
        return self.oniom_layer
    def is_link_atom_host(self):
        return self.link_atom_host

    def set_number(self,nr):
        self.number = nr
    def set_mulliken_charge(self,mq):
        self.mulliken_charge = mq
    def set_mulliken_spin_pop(self,ms):
        self.mulliken_spin_pop = ms
    def set_oniom_layer(self,ol):
        self.oniom_layer = ol
    def set_link_atom_host(self,lah):
        self.link_atom_host = lah

class Geometry(object):
    atoms = [] # a list of atom objects
    n_atoms = 0 # number of atoms
    scf_energy = 0.0 # SCF Energy in a.u.
    oniom_energy = None # ONIOM extrapolated energy in a.u.
    in_scan = False
    scan_point = None
    in_irc = False
    irc_point_number = None
    irc_net_reaction_coordinate = None
    irc_path_number = None
    geom_converged = False
    h_lah_atoms = [] # a list of atom objects from H-layer or link atom hosts

    # the class "constructor"
    def __init__(self, atoms):
        self.atoms = atoms
        self.n_atoms = len(self.atoms)
        
    def get_atoms(self):
        return self.atoms
    def get_n_atoms(self):
        return self.n_atoms
    def get_scf_energy(self):
        return self.scf_energy
    def get_oniom_energy(self):
        return self.oniom_energy
    def get_scan_point(self):
        return self.scan_point
    def get_irc_point_number(self):
        return self.irc_point_number
    def get_irc_net_reaction_coordinate(self):
        return self.irc_net_reaction_coordinate
    def get_irc_path_number(self):
        return self.irc_path_number
    def get_geom_converged(self):
        return self.geom_converged
    def get_h_lah_atoms(self):
        return Geometry(self.h_lah_atoms)
    
    def set_atoms(self,atoms):
        self.atoms = atoms
    def set_n_atoms(self):
        self.n_atoms = len(self.atoms)
    def set_scf_energy(self,scf):
        self.scf_energy = scf
    def set_oniom_energy(self,oniom_e):
        self.oniom_energy = oniom_e        
    def set_in_scan(self,in_scan):
        self.in_scan = in_scan
    def set_scan_point(self,pt):
        self.scan_point = pt
    def set_in_irc(self,in_irc):
        self.in_irc = in_irc
    def set_irc_point_number(self,irc_pt):
        self.irc_point_number = irc_pt
    def set_irc_net_reaction_coordinate(self,irc_net_coord):
        self.irc_net_reaction_coordinate = irc_net_coord
    def set_irc_path_number(self,irc_path_nr):
        self.irc_path_number = irc_path_nr
    def set_geom_converged(self,is_conv):
        self.geom_converged = is_conv
    def set_h_lah_atoms(self):
        h_lah_atoms = []
        for atom in self.atoms:
            if (atom.get_oniom_layer() == 'H') or atom.is_link_atom_host() :
                h_lah_atoms.append(atom)
        self.h_lah_atoms = h_lah_atoms
    def set_mulliken_charges(self,mull_q_list):
        if self.n_atoms == len(mull_q_list):
            for atom, q in zip(self.atoms, mull_q_list):
                atom.set_mulliken_charge(q)
        else:
            print("Error attempting to ascribe atom Mulliken charges \n")
            print("the geometry has: ", self.n_atoms, " atoms \n")
            print("while provided list of atom charges has length: ", len(mull_q_list))
            exit(1)
    def set_mulliken_spin_pops(self,mull_spop_list):
        if self.n_atoms == len(mull_spop_list):
            for atom, s in zip(self.atoms, mull_spop_list):
                atom.set_mulliken_spin_pop(s)
        else:
            print("Error attempting to ascribe atom Mulliken spin populations \n")
            print("the geometry has: ", self.n_atoms, " atoms \n")
            print("while provided list of atom spin pops has length: ", len(mull_spop_list))
            exit(1)
        
    def print_xyz(self):
        print( str(self.n_atoms) )
        if self.oniom_energy:
            head_line = 'oniom energy: ' + str(self.oniom_energy)
        else:
            head_line = 'scf energy: ' + str(self.scf_energy)
        if self.in_scan:     
            head_line = ' Scan point ' + str(self.scan_point) + ' ' + head_line
        elif self.in_irc:
            head_line = ' IRC net coordinate: ' + str(self.irc_net_reaction_coordinate) + ' ' + head_line            
        print(head_line)
        for atom in self.atoms:
            ele = atom.get_symbol()
            at_coord = atom.get_coords()
            line = ele + '\t' +\
                '{:06.6f}'.format(at_coord[0]) + '     ' +\
                '{:06.6f}'.format(at_coord[1]) + '     ' +\
                '{:06.6f}'.format(at_coord[2])
            print(line)       

    def print_xyzs(self):
        print( str(self.n_atoms) )
        if self.oniom_energy:
            head_line = 'oniom energy: ' + str(self.oniom_energy)
        else:
            head_line = 'scf energy: ' + str(self.scf_energy)
        if self.in_scan:     
            head_line = ' Scan point ' + str(self.scan_point) + ' ' + head_line
        elif self.in_irc:
            head_line = ' IRC net coordinate: ' + str(self.irc_net_reaction_coordinate) + ' ' + head_line            
        print(head_line)
        for atom in self.atoms:
            ele = atom.get_symbol()
            at_coord = atom.get_coords()
            spin_pop = atom.get_mulliken_spin_pop()
            line = ele + '\t' +\
                '{:06.6f}'.format(at_coord[0]) + '     ' +\
                '{:06.6f}'.format(at_coord[1]) + '     ' +\
                '{:06.6f}'.format(at_coord[2]) + '     ' +\
                '{:06.6f}'.format(spin_pop)    
            print(line)
            
        
def log_read_inp_geo(file, flag_line):
    """
    Reads initial geometry from the Gaussian output
    section marked with "Symbolic Z-matrix:" or "Symbolic Z-Matrix:"
    Parameters
    ----------
    file : log file (file object)
    flag_line : (string) 'Symbolic Z-matrix:' or 'Symbolic Z-Matrix:'
    Returns
    -------
    geometry (Geometry object) or string "EOF"
    """
    file.seek(0)
    atoms = []
    j = 0
    while True:
        a = file.readline()
        if not a:
            file.seek(0)
            return "EOF"
        match_flag=re.search(flag_line,a)
        if match_flag:
            while True:
                lah = False
                o_layer = None
                a = file.readline()
                a_list = a.split()
                a_len = len(a_list) 
                if a_len == 0:
                    break    
                elif a_list[0] == 'Charge':
                    pass
                else:
                    j += 1
                    at_symbol = a_list[0].split('-')[0]
                    at_number = symbol_at_num[at_symbol]
                    if (a_list[1] == '-1') or (a_list[1] == '0'):
                        x = eval(a_list[2])
                        y = eval(a_list[3])
                        z = eval(a_list[4])
                        if a_len > 5:
                            o_layer = a_list[5]
                        if a_len > 6:
                            lah = True
                    else:
                        x = eval(a_list[1])
                        y = eval(a_list[2])
                        z = eval(a_list[3])
                        if a_len > 4:
                            o_layer = a_list[4]
                        if a_len > 5:
                            lah = True
                    at_coords = [x, y, z]
                    atom = Atom(at_number, at_coords)
                    atom.set_oniom_layer(o_layer)
                    atom.set_link_atom_host(lah)
                    atom.set_number(j)
                    atoms.append(atom)
            geometry = Geometry(atoms)
            file.seek(0)
            return geometry


def log_read_geo(file, flag_line, n_atoms):
    """
    Reads atomic coordinates [A] in input orientation from the Gaussian output
    section marked with content of flag_line ("Input orientation:" or "Standard orientation:")
    Parameters
    ----------
    file : log file (file object)
    flag_line : line marking the geometry (string)
    n_atoms : number of atoms in the geometry to be read (int)
    Returns
    -------
    geometry (Geometry object) or string "EOF"
    """
    atoms = []
    fortran_line = ff.FortranRecordReader('(I7, I11, I16, 3F12.6)')
    while True:
        a = file.readline()
        if not a:
            return "EOF"
        match_flag=re.search(flag_line,a)
        if match_flag:            
            for i in range(4):
                file.readline()
            for i in range(n_atoms):
                a = file.readline()
                line = fortran_line.read(a)
                at_number = line[1]
                at_coords = [line[3], line[4], line[5]]
                atom = Atom(at_number, at_coords)
                atom.set_number(i+1)
                atoms.append(atom)
            geometry = Geometry(atoms)
            return geometry


def read_geo_scf_oniom_e(input_f, flag, n_atoms, ONIOM):
    """ from Gaussian output file read geometry and its SCF and ONIOM energies 
    INPUT: input_f - file object (Gaussian output file)
           flag - string, either "Input orientation:" or "Standard orientation:"
           n_atoms - number or atoms in the geoetry to read (int)
    OUTPUT: temp_geo - Geometry object or "EOF" string  """
    temp_geo = log_read_geo(input_f, flag, n_atoms)
    scf_e = None
    oniom_e = None
    if temp_geo != "EOF":
        scf_e = log_read_scf(input_f)
        temp_geo.set_scf_energy(scf_e)
        if ONIOM:
            oniom_e = log_read_oniom_e(input_f)
            temp_geo.set_oniom_energy(oniom_e)
    return temp_geo


def read_geo_mm_e(input_f, flag, n_atoms):
    """ from Gaussian output file read MM optimized geometry and its MM energy 
    INPUT: input_f - file object (Gaussian output file)
           flag - string, either "Input orientation:" or "Standard orientation:"
           n_atoms - number or atoms in the geoetry to read (int)
    OUTPUT: temp_geo - Geometry object or "EOF" string  """
    mm_e = log_read_mm_e(input_f)
    temp_geo = log_read_geo(input_f, flag, n_atoms)
    if temp_geo != "EOF":        
        temp_geo.set_scf_energy(mm_e)
        return temp_geo
    else:
        return "EOF" 


def log_read_mulliken(file):
    """
    read Mulliken charges and (if present) spin populations from Gaussian log file

    Parameters
    ----------
    file : - file object (Gaussian output file)
        
    Returns
    -------
    mull_q : list of floats
        a list of Mulliken atomic partial charges
    mull_s : list of floats
        a list of Mulliken atomic spin populations 

    """
    mull_q = []
    mull_s = []
    a = file.readline()
    a_list = a.split()
    a_len = len(a_list)
    if a_len == 2:
        s_present = True
    elif a_len == 1:
        s_present = False
    else:
        print("Error while trying to read Mulliken charges (and) spin populations \n")
        exit(1)
    while True:
        a = file.readline()
        if a[0:4] == " Sum":
            break    
        else:
            a_list = a.split()
            q = eval(a_list[2])
            mull_q.append(q)
            if s_present:
                s = eval(a_list[3])
                mull_s.append(s)
    return (mull_q, mull_s)


def step_line_tuple(x,y,z,w):
    step_tuple = namedtuple("step_tuple", "step_nr step_max scan_point scan_max")
    return step_tuple(*step_tuple(x,y,z,w))


def log_read_step_number_line(file):
    """
    From Gaussian log file reads a line "Step number   X out of a maximum of   Y on scan point     Z out of    W"
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    a tuple (step_nr, step_max, scan_point, scan_max)
    """
    # w pliku file wyszukaj odp linii
    flag_line = "Step number"
    while True:
        a = file.readline()
        if not a:
            break
        match_flag=re.search(flag_line,a)
        if match_flag:
            a_split = a.split() 
            step_nr = eval(a_split[2])
            step_max = eval(a_split[8])
            scan_point = eval(a_split[12])
            scan_max = eval(a_split[15])
            slt = step_line_tuple(step_nr,step_max,scan_point,scan_max)
            return slt
    return None


def irc_point_tuple(x,y,z):
    irc_tup = namedtuple("irc_tup", "path_nr point_nr net_reaction_coord")
    return irc_tup(*irc_tup(x,y,z))


def log_read_irc_data(file, last):
    """
    read IRC information from the Gaussian log file

    Parameters
    ----------
    file : file object
        gaussian log file
    last : int
        number of a last IRC point to be read from this file

    Returns
    -------
    irc_tup : namedtuple
        a tuple with: path_nr, point_nr, net_reaction_coord

    """
    flag_line = "Point Number:"
    while True:
        a = file.readline()
        if not a:
            break
        match_flag=re.search(flag_line,a)
        if match_flag:
            a_split = a.split() 
            point_nr = eval(a_split[2])
            path_nr = eval(a_split[5])
            file.readline()
            a = file.readline()
            a_split = a.split()
            net_reaction_coord = eval(a_split[8])
            irc_tup = irc_point_tuple(path_nr,point_nr,net_reaction_coord)
            if last == None:
                return irc_tup
            if last != None and point_nr <= int(last):
                return irc_tup
            else:
                break
    return None


def log_irc_or_scan(file):
    """
    In a Gaussian log file find if the calculations are IRC or scan
    If none of the above, returns "LAST"
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    a string "SCAN" or "IRC" or "LAST"
    """        
    file.seek(0)
    flag_irc = "IRC-IRC-IRC-IRC"
    flag_scan = "Step number"
    while True:
        a = file.readline()
        if not a:
            file.seek(0)
            return "LAST"
        match_irc = re.search(flag_irc,a)
        match_scan = re.search(flag_scan,a)
        if match_irc:
            file.seek(0)
            return "IRC"
        elif match_scan:
            file.seek(0)
            return "SCAN"


def log_is_ONIOM(file):
    """
    In a Gaussian log file find if the calculations are of ONIOM type
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    BOOL: True or False
    """        
    file.seek(0)
    flag_oniom = "ONIOM: generating point"
    flag_scf = "SCF Done:"
    while True:
        a = file.readline()
        if not a:
            file.seek(0)
            return False
        match_oniom = re.search(flag_oniom,a)
        match_scf = re.search(flag_scf,a)
        if match_oniom:
            file.seek(0)
            return True
        elif match_scf:
            file.seek(0)
            return False


def log_is_MM(file):
    """
    In a Gaussian log file find if the calculations are of MM type
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    BOOL: True or False
    """        
    file.seek(0)
    flag_amber = " AMBER calculation of energy"
    flag_oniom = "ONIOM: generating point"
    while True:
        a = file.readline()
        if not a:
            file.seek(0)
            return False
        match_amber = re.search(flag_amber,a)
        match_oniom = re.search(flag_oniom,a)
        if match_amber:
            file.seek(0)
            return True
        elif match_oniom:
            file.seek(0)
            return False
            
      
def log_read_scf(file):
    """
    In a Gaussian log file find SCF energy
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    a float: scf_energy or None
    """ 
    flag_no_conv = ">>>>>>>>>> Convergence criterion not met"
    flag_scf = "SCF Done:"
    while True:
        a = file.readline()
        if not a:
            return None
        match_no_conv = re.search(flag_no_conv,a)
        match_scf = re.search(flag_scf,a)
        if match_no_conv:
            file.readline() # skip the non-converged SCF Energy
        elif match_scf:
            a_split = a.split()
            scf_energy = eval( a_split[4] )
            return scf_energy


def log_read_oniom_e(file):
    """
    In a Gaussian log file find ONIOM energy
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    a float: oniom_energy
    """ 
    flag_oniom_e = "ONIOM: extrapolated energy ="
    while True:
        a = file.readline()
        if not a:
            break
        match_oniom_e = re.search(flag_oniom_e,a)
        if match_oniom_e:
            a_split = a.split()
            oniom_energy = eval( a_split[4] )
            return oniom_energy


def log_read_mm_e(file):
    """
    In a Gaussian log file find MM energy
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    a float: mm_energy or None
    """ 
    flag_e_f_class = "Energy per function class:"
    flag_e = "Energy="
    while True:
        a = file.readline()
        if not a:
            return None
        match_class = re.search(flag_e_f_class,a)
        if match_class:
            while True:
                a = file.readline()
                if not a:
                    return None
                match_e = re.search(flag_e,a)
                if match_e:
                    a_split = a.split()
                    mm_energy = eval( a_split[1] )
                    return mm_energy                    
                                        
        
def is_geom_converged(file):
    """
    In a Gaussian log file find if optimization completed or not
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    BOOL: True or False
    """ 
    flag_conv = "Optimization completed"
    flag_end_section = "GradGradGradGrad"
    while True:
        a = file.readline()
        if not a:
            break
        match_conv = re.search(flag_conv,a)
        if match_conv:
            return True
        else:
            match_end_section = re.search(flag_end_section,a)
            if match_end_section:
                return False


def is_irc_converged(file):
    """
    In a Gaussian log file find if IRC point optimization completed or not
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    BOOL: True or False
    """
    flag_conv = "Delta-x Convergence Met"
    flag_end_section = "Calculating another point on the path."
    while True:
        a = file.readline()
        if not a:
            break
        match_conv = re.search(flag_conv,a)
        match_end_section = re.search(flag_end_section,a)
        if match_end_section:
            return False
        elif match_conv:
            return True    


def scan_file(file, text_flag):
    """
    scan the whole file of occurance of text_flag

    Parameters
    ----------
    file : file objects
        text file
    text_flag : string
        a pattern to find

    Returns
    -------
    positions : list
        positions (just after) of the pattern found in the file 

    """
    file.seek(0)
    positions = []
    while True:
        a = file.readline()
        if not a:
            break
        match = re.search(text_flag,a)
        if match:
            positions.append(file.tell())
    file.seek(0)
    return positions


def find_first_infile(file, text_flag):
    """
    finds first instance of text_flag in the file

    Parameters
    ----------
    file : file objects
        text file
    text_flag : string
        pattern to find

    Returns
    -------
    positions : list
        empty list if pattern not found
        1-element list with position, if pattern found

    """
    file.seek(0)
    positions = []
    while True:
        a = file.readline()
        if not a:
            break
        match = re.search(text_flag,a)
        if match:
            positions.append(file.tell())
            file.seek(0)
            return positions        


def print_help():  
    help_text = """
This script reads Gaussian log file and outputs xyz geometry/ies.

For geometry 1-D scan or IRC calculations optimized (or last) geometries 
for points along the scanned coordinate or IRC are output. A png file with a plot
of ONIOM or SCF energies along the profile is generated.

The script expects 2 or 3 arguments: 
    #1 log-file-name or, for irc_f, name of a file specifying filenames and directions of irc calc.
    #2 type of extraction - one from amoung: scan, irc, irc_f, all, last, nr, lasts, nrs
    #3 (for IRC) file name with SP/FREQ calculations for the TS from which IRC calculations started
       or for NR 1-based number of the structure to be extracted
    """
    
    print(help_text) 