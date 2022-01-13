#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
a collection of classes and functions for gaussian2xyz

Authors: Tomasz Borowski
"""

import numpy as np
import re
from collections import namedtuple

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


class Atom(object):
    atomic_number = 0
    symbol = ""
    number = 0
    coordinates = np.empty((1, 3), dtype=np.float64)  # coordinates in Angstroms
    mulliken_charge = None
    mulliken_spin_pop = None

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

    def set_number(self,nr):
        self.number = nr
    def set_mulliken_charge(self,mq):
        self.mulliken_charge = mq
    def set_mulliken_spin_pop(self,ms):
        self.mulliken_spin_pop = ms    


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
       
def log_read_geo(file):
    """
    Reads atomic coordinates [A] from the Gaussian output
    section marked with "Input orientation:" 
    Parameters
    ----------
    file : log file (file object)
    Returns
    -------
    geometry (Geometry object) or string "EOF"
    """
    # w pliku file wyszukaj odp linii
    flag_line = "Input orientation:"
    atoms = []
    j=0
    while True:
        a = file.readline()
        if not a:
            return "EOF"
            print(j)
        match_flag=re.search(flag_line,a)
        if match_flag:
            j=j+1
            for i in range(4):
                file.readline()
            while True:
                a = file.readline()
                if a[1] == "-":
                    break    
                else:
                    a_split = a.split()
                    at_number = eval( a_split[1] )
                    at_coords = [eval(a_split[3]), eval(a_split[4]), eval(a_split[5])]
                    atom = Atom(at_number, at_coords)
                    atoms.append(atom)
            geometry = Geometry(atoms)
            return geometry

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
#            return (step_nr, step_max, scan_point, scan_max)
            slt = step_line_tuple(step_nr,step_max,scan_point,scan_max)
            return slt
    return None


def irc_point_tuple(x,y,z):
    irc_tup = namedtuple("irc_tup", "path_nr point_nr net_reaction_coord")
    return irc_tup(*irc_tup(x,y,z))


def log_read_irc_data(file, last):
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
        match_end_section = re.search(flag_end_section,a)
        if match_end_section:
            return False
        elif match_conv:
            return True

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
    #flag_end_section = "IRC-IRC-IRC-IRC"
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
