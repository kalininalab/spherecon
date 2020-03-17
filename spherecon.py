#!/usr/bin/python3
import os
import sys
import getopt
import math
import numpy
import gzip
from Bio.Data import SCOPData as scd

threeToOne = scd.protein_letters_3to1

oneToThree = {'C':'CYS', 
'D':'ASP',
'S':'SER',
'V':'VAL',
'Q':'GLN',
'K':'LYS',
'P':'PRO',
'T':'THR',
'F':'PHE',
'A':'ALA',
'H':'HIS',
'G':'GLY',
'I':'ILE',
'E':'GLU',
'L':'LEU',
'R':'ARG',
'W':'TRP',
'N':'ASN',
'Y':'TYR',
'M':'MET',
'X':'UNK'}

vdw_radius = {"C":1.7,"O":1.52,"N":1.55,"F":1.47,"P":1.8,"S":1.8,"X":1.7,'A':1.85,'I':1.98,'B':1.85,'R':2.0, 'V':1.53}

radii_map = {'CYS': 3.029584442714009, 'SEC': 3.029584442714009,'ILE': 3.3236043438778844, 'SER': 2.942862799880502, 'GLN': 3.392265173640708, 'LYS': 3.4323769248759275, 'TRP': 4.020778501571494, 'PRO': 3.1681967871564782, 'THR': 3.1209640826827605, 'PHE': 3.7193696421024685, 'ALA': 2.8009640812798096, 'GLY': 2.5743877263600927, 'HIS': 3.534660601593296, 'ASN': 3.2435252870320834, 'LEU': 3.3236043438778844, 'ARG': 3.631344304999469, 'ASP': 3.236792124221129, 'VAL': 3.1681967871564782, 'GLU': 3.3861111330098956, 'TYR': 3.8021338661518342, 'MET': 3.351107764635289}

avg_cent_side = {'ILE': (-0.90613294951389123, -1.2696649328274769, 0.75787678394555402), 'GLN': (-1.2338776755913703, -1.6271318863182545, 0.98235396198659397), 'GLY': (-0.0047299661865784631, -0.0088845691075302401, 9.8567539791531887e-05), 'GLU': (-1.2149070583116599, -1.5273224890561412, 1.0878295037481032), 'CYS': (-0.58968190337597537, -0.89958163307781258, 0.62776825003201231), 'HIS': (-0.9776890826599115, -1.6575223413455686, 0.78038430859923047), 'SER': (-0.59170713881015236, -0.65044117400421442, 0.79631448621746626), 'LYS': (-1.4201311703581181, -1.8190188812624906, 1.0949225828431293), 'PRO': (0.21186055970214385, -0.81932530473289122, 1.0393498184079177), 'SEC': (-0.83941591976390406, -0.71433578942505138, 0.85735910033472262), 'ASN': (-0.84857508137952831, -1.3074529598975138, 0.7776964830504921), 'VAL': (-0.79991957561390237, -0.94684856115170413, 0.62104228622222646), 'THR': (-0.65195583478889552, -0.90724937534582129, 0.82850637927897564), 'ASP': (-0.90583763237409776, -1.1848068294814498, 0.82225350991223289), 'TRP': (-1.528171801920162, -1.7410619994461729, 0.89554186538385805), 'UNK': (-0.3674831972563582, -0.32646621880413212, 0.49726141188550765), 'PHE': (-1.1835851942281017, -1.7644149833489573, 0.692845774001975), 'ALA': (-0.40699060927186675, -0.36885898636123338, 0.49092580988001461), 'MET': (-1.1860699223797722, -1.5396420763646954, 0.84853591931547168), 'LEU': (-1.0557131262899473, -1.4510720809764079, 0.67350078632184163), 'ARG': (-1.7179619489063336, -2.0937800414713896, 1.3885442733060744), 'TYR': (-1.2761868751958905, -1.9287097734708729, 0.69640859762760676)}

unknown_avg_cent_side = (-0.866137516594, -1.16682367833, 0.769774990962) #average vector over all residue types
unknown_thresh = 8.0
unknown_angle_thresh = -0.8
unknown_rad = 3.23342916549

avg_cent = {'ILE': (-0.51939923132192789, -0.42924346469309471, 0.57478473070875846), 'GLN': (-0.80834986146513732, -0.76972077063003586, 0.68335543707029034), 'GLY': (0.003710108790125123, 0.74850253699086222, -0.011107678632392972), 'GLU': (-0.78586343234711975, -0.68596611605471958, 0.72074313020193725), 'CYS': (-0.25453284884232952, 0.025290152329503414, 0.46347585732429047), 'HIS': (-0.67516717655691016, -0.87046086894289898, 0.61423651641973831), 'SER': (-0.25186295200626596, 0.15435437507216301, 0.50821742795023583), 'LYS': (-0.93157547397808582, -0.90270366121918721, 0.7564363100362167), 'PRO': (0.19092956143389894, -0.054932017506649618, 0.7190732935459736), 'SEC': (-0.4041486429034104, 0.19276231247130715, 0.31022206884216508), 'ASN': (-0.52337772135489558, -0.46464652365549852, 0.5839222978918297), 'VAL': (-0.40625381747595429, -0.12555581649592371, 0.49308718269365215), 'THR': (-0.3458279994581282, -0.10156321384665247, 0.56845475732264672), 'ASP': (-0.5395255947207902, -0.37285940531331668, 0.57492539764983719), 'TRP': (-1.1801359800054108, -1.1546490468404469, 0.74010431975141766), 'UNK': (-0.16216122024915333, 0.48606836524146568, 0.13849611887759408), 'PHE': (-0.83581175433622168, -1.0098420145081217, 0.57161519965229146), 'ALA': (-0.11519817930399062, 0.43429992078978386, 0.21401964289355344), 'MET': (-0.70928160313770094, -0.59184549532823161, 0.58483387893500882), 'LEU': (-0.62691092485685773, -0.54401570724965409, 0.4676475183702658), 'ARG': (-1.2457704954369544, -1.2737121830387164, 1.037835446531018), 'TYR': (-0.9381954559286525, -1.1994693553734199, 0.57799294311076543)}

threshs_alpha ={
'CYS' : 7.25, 
'ASP' : 7.25,
'SER' : 7.0,
'VAL' : 7.5,
'GLN' : 7.75,
'LYS' : 8.0,
'PRO' : 7.5,
'THR' : 7.0,
'PHE' : 8.25,
'ALA' : 7.0,
'HIS' : 8.0,
'GLY' : 6.75,
'ILE' : 7.75,
'GLU' : 7.5,
'LEU' : 7.75,
'ARG' : 8.5,
'TRP' : 9.0,
'ASN' : 7.25,
'TYR' : 8.75,
'MET' : 8.0
}

angle_threshs_alpha = {
'CYS' : -0.85,
'ASP' : -0.90,
'SER' : -0.85,
'VAL' : -0.90,
'GLN' : -1.0,
'LYS' : -1.0,
'PRO' : -1.0,
'THR' : -1.0,
'PHE' : -0.90,
'ALA' : -0.85,
'HIS' : -0.95,
'GLY' : -1.0,
'ILE' : -0.85,
'GLU' : -1.0,
'LEU' : -1.0,
'ARG' : -1.0,
'TRP' : -0.6,
'ASN' : -1.0,
'TYR' : -0.85,
'MET' : -1.0
}

#learned values from golden train set
threshs = {'CYS' : 7.0, 'SEC': 7.0, 'GLN' : 7.25, 'ASP' : 6.75, 'ASX': 6.75, 'SER' : 6.5, 'VAL' : 7.5, 'LYS' : 7.5, 'ASN' : 6.75, 'PRO' : 7.0, 'THR' : 6.75, 'PHE' : 7.75, 'ALA' : 7.0, 'HIS' : 7.25, 'GLY' : 6.75, 'ILE' : 7.75 , 'LEU' : 7.75, 'ARG' : 7.75, 'TRP' : 8.0, 'GLU' : 7.0, 'GLX': 7.0, 'TYR' : 7.75, 'MET' : 7.5, 'UNK': unknown_thresh}

angle_threshs = {'CYS' : -0.85, 'SEC' : -0.85, 'GLN' : -0.95, 'ASP' : -0.95, 'ASX' : -0.95, 'SER' : -0.9, 'VAL' : -0.85, 'LYS' : -0.9, 'ASN' : -0.95, 'PRO' : -0.85, 'THR' : -0.85, 'PHE' : -0.85, 'ALA' : -0.85, 'HIS' : -0.8, 'GLY' : -1.0, 'ILE' : -0.9, 'LEU' : -0.9, 'ARG' : -0.95, 'TRP' : -0.95, 'GLU' : -0.95, 'GLX' : -0.95, 'TYR' : -0.85, 'MET' : -1.0, 'UNK': unknown_angle_thresh}

dist_matrix_threshs = {'None':{'GLY': 8.0, 'LYS': 8.5, 'TRP': 8.0, 'ARG': 8.0, 'GLU': 8.0, 'GLN': 8.0, 'LEU': 8.0, 'PHE': 8.0, 'ILE': 8.0, 'VAL': 8.0, 'HIS': 8.0, 'THR': 8.0, 'ASN': 8.0, 'SER': 8.0, 'MET': 9.0, 'ALA': 8.0, 'TYR': 8.0, 'ASP': 8.0, 'CYS': 8.0, 'PRO': 8.0},
                        0.0:{'GLY': 7.5, 'LYS': 8.0, 'TRP': 7.5, 'ARG': 7.5, 'GLU': 7.5, 'GLN': 7.5, 'LEU': 7.5, 'PHE': 7.5, 'ILE': 7.5, 'VAL': 7.5, 'HIS': 7.5, 'THR': 7.5, 'ASN': 7.5, 'SER': 7.5, 'MET': 8.0, 'ALA': 7.5, 'ASP': 7.5, 'PRO': 7.5, 'TYR': 7.5, 'CYS': 7.5},
                        0.1:{'GLU': 7.75, 'TYR': 7.75, 'PHE': 7.5, 'THR': 7.75, 'LEU': 7.5, 'GLN': 7.75, 'ILE': 7.5, 'ARG': 7.75, 'GLY': 7.5, 'LYS': 8.0, 'ASN': 7.75, 'ALA': 7.5, 'ASP': 7.75, 'SER': 7.5, 'PRO': 7.75, 'VAL': 7.5, 'CYS': 10.0, 'HIS': 7.75, 'MET': 7.5, 'TRP': 9.0},
                        0.2:{'MET': 9.25, 'GLY': 8.75, 'ILE': 8.5, 'LEU': 8.5, 'GLN': 8.75, 'PRO': 8.75, 'ASP': 8.75, 'ALA': 8.5, 'ASN': 8.75, 'THR': 8.75, 'VAL': 8.5, 'ARG': 8.75, 'TRP': 8.5, 'SER': 8.5, 'HIS': 8.75, 'PHE': 8.5, 'TYR': 9.0, 'GLU': 9.5, 'LYS': 9.0, 'CYS': 9.0},
                        0.3:{'MET': 9.25, 'GLY': 8.5, 'ILE': 8.5, 'LEU': 8.5, 'GLN': 8.75, 'PRO': 8.75, 'ASP': 8.75, 'ALA': 8.5, 'ASN': 8.75, 'THR': 8.75, 'VAL': 8.5, 'ARG': 8.75, 'TRP': 9.0, 'SER': 8.5, 'HIS': 9.0, 'PHE': 8.5, 'TYR': 9.5, 'GLU': 8.75, 'LYS': 9.25, 'CYS': 8.25},
                        0.4:{'MET': 9.75, 'GLY': 8.25, 'ILE': 8.25, 'LEU': 8.0, 'GLN': 8.5, 'PRO': 8.25, 'ASP': 8.25, 'ALA': 8.25, 'ASN': 8.25, 'THR': 8.25, 'VAL': 9.0, 'ARG': 8.25, 'TRP': 9.25, 'SER': 8.25, 'HIS': 8.25, 'PHE': 8.0, 'TYR': 8.25, 'GLU': 8.5, 'LYS': 9.0, 'CYS': 8.25},
                        0.5:{'MET': 9.25, 'GLY': 8.75, 'ILE': 8.5, 'LEU': 8.5, 'GLN': 9.0, 'PRO': 8.75, 'ASP': 8.75, 'ALA': 8.5, 'ASN': 9.0, 'THR': 9.0, 'VAL': 8.5, 'ARG': 8.75, 'TRP': 9.0, 'SER': 8.5, 'HIS': 8.75, 'PHE': 9.0, 'TYR': 8.75, 'GLU': 9.0, 'LYS': 9.25, 'CYS': 8.5},
                        0.6:{'MET': 10.25, 'LYS': 9.25, 'LEU': 8.5, 'ASP': 9.0, 'GLU': 9.5, 'ILE': 8.5, 'ALA': 8.5, 'ARG': 9.0, 'GLY': 8.75, 'VAL': 8.5, 'SER': 8.75, 'THR': 8.75, 'TYR': 9.0, 'ASN': 9.0, 'GLN': 9.25, 'HIS': 8.75, 'PRO': 9.0, 'PHE': 9.0, 'CYS': 10.0, 'TRP': 9.5},
                        0.7:{'VAL': 9.75, 'THR': 9.75, 'LYS': 10.5, 'CYS': 9.5, 'ASN': 10.0, 'ILE': 10.0, 'SER': 9.75, 'MET': 9.5, 'PRO': 10.0, 'ALA': 9.5, 'LEU': 9.75, 'HIS': 10.0, 'TYR': 9.75, 'GLN': 10.0, 'GLY': 9.75, 'ARG': 10.25, 'GLU': 10.25, 'PHE': 9.5, 'ASP': 10.25, 'TRP': 9.75},
                        0.8:{'ASN': 10.0, 'SER': 9.5, 'ILE': 9.25, 'VAL': 9.25, 'PRO': 10.0, 'ARG': 9.75, 'GLN': 9.75, 'GLY': 9.5, 'LEU': 10.0, 'LYS': 10.5, 'PHE': 9.5, 'TRP': 9.5, 'GLU': 10.25, 'ASP': 10.25, 'TYR': 9.5, 'THR': 9.75, 'CYS': 9.5, 'ALA': 9.5, 'HIS': 10.0, 'MET': 9.5},
                        0.9:{'MET': 9.25, 'SER': 8.75, 'ILE': 9.0, 'ALA': 8.75, 'GLU': 10.0, 'ARG': 9.25, 'GLN': 9.0, 'LEU': 8.5, 'PHE': 8.5, 'PRO': 9.25, 'VAL': 8.5, 'ASN': 9.25, 'THR': 8.75, 'ASP': 9.25, 'HIS': 9.0, 'GLY': 8.75, 'TYR': 8.75, 'TRP': 8.5, 'LYS': 9.5, 'CYS': 8.75}
}


def calcVol(r,cos):
    vol = (2.0/3.0)*math.pi*(r**3.0)*(1.0-cos)
    return vol

def sissCorrectionVol(siss,r,cos,target_rad):
    v = calcVol(r,cos) - calcVol(target_rad,cos)
    #print(r,cos,target_rad,v)
    return (v-siss)/v

def parsePDB(input_file,chains,c_alpha,target_atom='CA',page=None):
    """
    Parses a PDB-file and takes all atomic coordinates of a specified chain.

    Input:
    input_file: String ; Path to a PDB file
    chain: String or None ; Chain identifier, if None is given, the first Chain found in the file is taken
    c_alpha: Boolean ; If True, then only C alpha atoms are taken

    Output:
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.   
    """
    if page == None:
        try:
            if input_file[-4:] == '.pdb' or input_file[-5:-1] == '.pdb':
                f = open(input_file,'rb')
                lines = f.readlines()
                f.close()
            elif input_file[-3:] == '.gz':
                f =gzip.open(input_file,'rb')
                lines = f.readlines()
                f.close()
            else:
                raise NameError('Invalid input file: %s' % input_file)
        except:
            raise NameError('Invalid input file: %s' % input_file)
    else:
        lines = page.split('\n')

    coordinate_map = {}

    x_total = 0.0
    y_total = 0.0
    z_total = 0.0
    n = 0.0

    if len(lines) == 0:
        print(input_file)

    for line in lines:
        try:
            line = line.decode('ascii')
        except:
            pass
        if len(line) > 5:        
            record_name = line[0:6].replace(" ","")
            if record_name == "ENDMDL":
                break
        #ignore short lines
        if len(line) > 20:
            atom_nr = line[6:11].replace(" ","")
            atom_name = line[12:16].replace(" ","")
            res_name = line[17:20].replace(" ","")
            
            #len(res_name) == 3 is only true for amino acid chains
            if len(res_name) == 3:
                if len(line) > 21:
                    chain_id = line[21]
                    res_nr = line[22:27].replace(" ","")
                
                #consider only lines with record name ATOM
                if record_name == "ATOM" or record_name == 'HETATM':
                    if record_name == 'HETATM':
                        if res_name not in threeToOne:
                            continue
                        else:
                            res_name = OneToThree[threeToOne[res_name][0]]
                    #if chain not given, take first one found
                    if chains == None:
                        chains = [chain_id]
                    if len(line) > 50:
                        #print(chain_id,chain)
                        #consider only lines with the correct chain id
                        if chain_id in chains:
                            if not chain_id in coordinate_map:
                                coordinate_map[chain_id] = {}
                            #if c_alpha is true, then take only C alpha atoms
                            if res_name != 'UNK':
                                if (not c_alpha) or atom_name == target_atom:
                                    if atom_name[0] != 'H' and atom_name[0] != 'D':
                                        x = float(line[30:38].replace(" ",""))
                                        y = float(line[38:46].replace(" ",""))
                                        z = float(line[46:54].replace(" ",""))
                                        if res_nr not in coordinate_map[chain_id]:
                                            coordinate_map[chain_id][res_nr]=[res_name,{}]
                                        coordinate_map[chain_id][res_nr][1][atom_nr] = (atom_name,x,y,z)
                                        x_total += x
                                        y_total += y
                                        z_total += z
                                        n += 1.0
    if n > 0.0:                                        
        protein_centroid = numpy.array([x_total/n,y_total/n,z_total/n])
    else:
        protein_centroid = numpy.array([0.0,0.0,0.0])
    #print(input_file)
    if len(coordinate_map) == 0:
        print('Warning: Chain identifier not in given PDB structure: ',input_file,chains)

    return coordinate_map,protein_centroid

def calcCentroidMap(coordinate_map,target_residues,c_alpha,double_unknown_mode = False,target_atom = 'CA'):
    centroid_map = {}
    if c_alpha:
        c1 = None
        c2 = None
        res_2 = None
        res_1 = None
        res_name_1 = None
        res_name_2 = None
        res_nr_1 = None
        res_nr_2 = None
        for chain in coordinate_map:
            centroid_map[chain] = {}
            for res in coordinate_map[chain]:
                try:
                    res_nr = int(res)
                except:
                    res_nr = int(res[:-1])
                c0 = c1
                c1 = c2
                res_nr_0 = res_nr_1
                res_nr_1 = res_nr_2
                res_nr_2 = res_nr
                res_0 = res_1
                res_1 = res_2
                res_2 = res
                res_name_1 = res_name_2
                res_name_2 = coordinate_map[chain][res_2][0]
                atomlist = coordinate_map[chain][res_2][1]
                (atomname,x,y,z) = list(atomlist.values())[0]
                if atomname != target_atom:
                    raise NameError('Only atom is not C alpha')
                c2 = numpy.array([x,y,z])           
                if res_0 != None and res_1 != None:
                    if res_nr_2 - res_nr_1 == 1 and res_nr_1 - res_nr_0 == 1:
                        avg_cent_side_vec = avg_cent_side[res_name_1]
                        avg_cent_vec = avg_cent[res_name_1]
                        if double_unknown_mode: #in double unknown mode, the residue types are not needed!
                            avg_cent_side_vec = unknown_avg_cent_side
                        side_centroid = predict_centroid(c0,c1,c2,avg_cent_side_vec)
                        centroid = predict_centroid(c0,c1,c2,avg_cent_vec)
                        
                        centroid_map[chain][res_1] = (side_centroid,centroid)

                    else:
                        centroid_map[chain][res_1] = (c1,c1)
                elif res_1 != None:
                    centroid_map[chain][res_1] = (c1,c1)
            if len(coordinate_map[chain]) > 0:
                centroid_map[chain][res_2] = (c2,c2)
    else:
        for chain in target_residues:
            centroid_map[chain] = {}
            for res in target_residues[chain]:
                if not res in coordinate_map[chain]:
                    print('Warning: given residue %s not found' % res)
                    continue
                atomlist = coordinate_map[chain][res][1]
                centroid = getCentroid(atomlist)
                centroid_map[chain][res] = centroid
    return centroid_map

def parseDM(input_file,add_neighbors=False):
    f = open(input_file,'r')
    lines = f.read().split('\n')
    f.close()

    dist_matrix = {}

    seq = ''

    for line in lines:
        if line == '':
            continue
        words = line.split()

        if words[0] == 'SEQ':
            seq += words[1]


        if ord(words[0][0]) < 18 or ord(words[0][0]) > 57:
            continue

        if len(words) < 3:
            continue

        res_1 = words[0]
        res_2 = words[1]
        d = float(words[2])

        dist_matrix[('A',res_1,'A',res_2)] = d
        dist_matrix[('A',res_2,'A',res_1)] = d

    res_name_map = {'A':{}}

    for pos,letter in enumerate(seq):
        res_name_map['A'][str(pos+1)] = oneToThree[letter]

    for res in res_name_map['A']:
        res2 = str(int(res)+1)
        if not res2 in res_name_map['A']:
            continue
        if add_neighbors:
            if not ('A',res,'A',res2) in dist_matrix:
                dist_matrix[('A',res,'A',res2)] = 6.92
                dist_matrix[('A',res2,'A',res)] = 6.92
    N = 0
    M = 0
    for pos1,res1 in enumerate(seq):
        res1 = str(pos1+1)
        for pos2,res2 in enumerate(seq):
            res2 = str(pos2+1)
            if abs(pos2-pos1) < 3:
                continue
            N+=1
            if (('A',res1,'A',res2)) in dist_matrix:
                M += 1

    completeness = M/N
    sparsity = 1.-completeness
    sparsity = min(0.9,round(sparsity,1))
    return dist_matrix,res_name_map,seq,sparsity

def calcDistMatrix(coordinate_map,centroid_map,target_residues,c_alpha,atom_mode = False):
    dist_matrix = {}
    if c_alpha:
        for chain in centroid_map:
            for res_1 in centroid_map[chain]:
                centroid_1 = centroid_map[chain][res_1][0]
                for chain2 in centroid_map:
                    for res_2 in centroid_map[chain2]:
                        if not (chain,res_1,chain2,res_2) in dist_matrix:
                            centroid_2 = centroid_map[chain2][res_2][0]
                            diff = centroid_2 - centroid_1
                            d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                            dist_matrix[(chain,res_1,chain2,res_2)] = d
                            dist_matrix[(chain2,res_2,chain,res_1)] = d
    elif not atom_mode:
        for chain in target_residues:
            for res_1 in target_residues[chain]:
                if not res_1 in centroid_map[chain]:
                    continue
                centroid = centroid_map[chain][res_1]
                res_name = coordinate_map[chain][res_1][0]
                thresh = threshs[res_name]
                dist_matrix[(chain,res_1)] = {}
                for chain2 in target_residues:
                    for res_2 in coordinate_map[chain2]:
                        dist_matrix[(chain,res_1)][(chain2,res_2)] = {}
                        atomlist = coordinate_map[chain2][res_2][1]
                        for atom in atomlist:
                            (atomname,x,y,z) = atomlist[atom]
                            diff = centroid - numpy.array([x,y,z])
                            d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                            #If one atom is further than 10A+threshhold, then ignore the whole residue
                            if d - 10.0 > thresh:
                                break

                            dist_matrix[(chain,res_1)][(chain2,res_2)][atom] = (d,atomname,x,y,z)
    else:
        for chain in target_residues:
            for res_1 in target_residues[chain]:

                res_name = coordinate_map[chain][res_1][0]
                atomlist_1 = coordinate_map[chain][res_1][1]
                thresh = threshs[res_name]
                dist_matrix[(chain,res_1)] = {}
                for chain2 in target_residues:
                    for res_2 in coordinate_map[chain2]:
                        dist_matrix[(chain,res_1)][(chain2,res_2)] = {}
                        atomlist_2 = coordinate_map[chain2][res_2][1]
                        for atom in atomlist_1:
                            dist_matrix[(chain,res_1)][(chain2,res_2)][atom] = {}
                            for atom2 in atomlist_2:
                                (atomname1,x1,y1,z1) = atomlist_1[atom]
                                (atomname2,x2,y2,z2) = atomlist_2[atom2]
                                diff = numpy.array([x1,y1,z1]) - numpy.array([x2,y2,z2])
                                d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)

                                dist_matrix[(chain,res_1)][(chain2,res_2)][atom][atom2] = (d,atomname2,x2,y2,z2)
    return dist_matrix

def sphere_intersection(R,r,d):
    if R < r:
        a = r
        r = R
        R = a
    if r+R <= d:
        return 0.0
    if d+r <= R:
        return (4.0/3.0)*math.pi*r**3.0
    sum1 = (R+r-d)**2
    sum2 = (d**2.0+2.0*d*r-3.0*r**2.0+2.0*d*R+6.0*r*R-3.0*R**2.0)
    si = math.pi*sum1*sum2/(12.0*d)
    #print('si',R,r,d,si)
    return si

def createRotMatrix(axis,cos):
    sin = (1.0-cos**2.0)**(1.0/2.0)
    if axis == "x":
        Rot = [[1.0,0.0,0.0],[0.0,cos,-sin],[0.0,sin,cos]]
    if axis == "y":
        Rot = [[cos,0.0,sin],[0.0,1.0,0.0],[-sin,0.0,cos]]
    if axis == "z":
        Rot = [[cos,-sin,0.0],[sin,cos,0.0],[0.0,0.0,1.0]]
    return numpy.matrix(Rot)

def createPlaneRotMatrix(plane,vec):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    if plane == 'xy':
        if z == 0.0:
            cos = 1.0
        elif y == 0.0:
            cos = 0.0
        else:
            cos = (((z**2.0/y**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('x',cos)
        if (vec*rot).A1[2]**2 > 0.0000001:
            rot = rot.T
    if plane == 'xz':
        if y == 0.0:
            cos = 1.0
        elif x == 0.0:
            cos = 0.0
        else:
            cos = (((y**2.0/x**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('z',cos)
        if (vec*rot).A1[1]**2 > 0.0000001:
            rot = rot.T
    if plane == 'yz':
        if x == 0.0:
            cos = 1.0
        elif z == 0.0:
            cos = 0.0
        else:
            cos = (((x**2.0/z**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('y',cos)
        if (vec*rot).A1[0]**2 > 0.0000001:
            rot = rot.T
    return rot

def createAxisRotMatrix(axis,vec):
    if axis == 'x':
        rot1 = createPlaneRotMatrix('xz',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('y',getCos('x',ivec))
        if (ivec*rot2).A1[2]**2 > 0.0000001:
            rot2 = rot2.T
    if axis == 'y':
        rot1 = createPlaneRotMatrix('xy',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('z',getCos('y',ivec))
        if (ivec*rot2).A1[0]**2 > 0.0000001:
            rot2 = rot2.T
    if axis == 'z':
        rot1 = createPlaneRotMatrix('yz',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('x',getCos('z',ivec))
        if (ivec*rot2).A1[1]**2 > 0.0000001:
            rot2 = rot2.T
    return rot1*rot2

def createRotAxisMatrix(axis,cos):
    sin = (1.0-cos**2.0)**(1.0/2.0)
    axis = axis/numpy.linalg.norm(axis)
    x = axis[0]
    y = axis[1]
    z = axis[2]
    m1 = [cos+x**2.0*(1.0-cos),x*y*(1.0-cos)-z*sin,x*z*(1.0-cos)+y*sin]
    m2 = [y*x*(1.0-cos)+z*sin,cos+y**2.0*(1.0-cos),y*z*(1.0-cos)-x*sin]
    m3 = [z*x*(1.0-cos)-y*sin,z*y*(1.0-cos)+x*sin,cos+z**2.0*(1.0-cos)]
    Rot = [m1,m2,m3]
    return numpy.matrix(Rot)

def gly_vector(n_v,c_v,ca_v):
    n_v = n_v - ca_v 
    c_v = c_v - ca_v
    rot = createRotAxisMatrix(c_v,-0.5)
    vec = (n_v*rot).A1
    return vec
    
def getCosAngle(vec1,vec2):
    n1 = numpy.linalg.norm(vec1)
    n2 = numpy.linalg.norm(vec2)
    norm = n1*n2
    dot = numpy.dot(vec1,vec2)
    if norm != 0.0:
        c = dot / norm
    else:
        return None
    # Take care of roundoff errors 
    c = min(c, 1.0) 
    c = max(-1.0, c)
    return c

def getCos(axis,vec):
    if axis == "x":
        other = numpy.array([1.0,0.0,0.0])
    if axis == "y":
        other = numpy.array([0.0,1.0,0.0])
    if axis == "z":
        other = numpy.array([0.0,0.0,1.0])
    n1 = numpy.linalg.norm(vec)
    n2 = numpy.linalg.norm(other)
    c = ((numpy.dot(vec,other)) / (n1 * n2))
    # Take care of roundoff errors 
    c = min(c, 1.0) 
    c = max(-1.0, c)
    return c

def nullTest(vec):
    if vec[0] == 0.0 and vec[1] == 0.0 and vec[2] == 0.0:
        return False
    return True

def predict_centroid(c0,c1,c2,avg_cent_vec):
    if nullTest(c0-c1) and nullTest(c2-c1):
        A = (c0-c1)/numpy.linalg.norm(c0-c1)
        B = (c2-c1)/numpy.linalg.norm(c2-c1)
        rx = createAxisRotMatrix("x",A)
        B_prime = (B*rx).A1
        r2 = createPlaneRotMatrix('xy',B_prime)
        ROT = rx*r2
        support = (B_prime*r2).A1
        if support[1] < 0.0:
            flip = createRotMatrix('x',-1.0)
            ROT = ROT*flip
        return ((avg_cent_vec*ROT.T).A1)+c1
    else:
        return c1

def getCentroid(atomlist):
    n = 0.0
    t_x = 0.0
    t_y = 0.0
    t_z = 0.0
    c_a = None
    c_b = None

    for atom in atomlist:
        (atomname,x,y,z) = atomlist[atom]
        if len(atomname) > 1:
            t_x += x
            t_y += y
            t_z += z
            n += 1.0
            
    if n > 0.0:
        centroid = numpy.array([t_x/n,t_y/n,t_z/n])
        
    else:
        for atom in atomlist:
            (atomname,x,y,z) = atomlist[atom]        
            t_x += x
            t_y += y
            t_z += z
            n += 1.0
            
        if n > 0.0:
            centroid = numpy.array([t_x/n,t_y/n,t_z/n])
        else:
            centroid = None
    return centroid

def produceOutput(siss_map,coordinate_map,output_file,res_name_map,atom_mode=False):
    """
    Input:
    siss_map: {String:float} ; 
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.
    """ 
    if len(siss_map) == 0:
        print('Warning: no output was produced')
        return

    lines = []    
    for chain in siss_map:
        for res in siss_map[chain]:
            if coordinate_map[chain] != None:
                res_name = coordinate_map[chain][res][0]
            else:
                res_name = res_name_map[chain][res]
            if not atom_mode:
                siss = siss_map[chain][res]
                lines.append("%s %s %s\t%s\n" % (res,chain,res_name,str(siss)))
            else:
                for atom in siss_map[chain][res]:
                    
                    siss = siss_map[chain][res][atom]
                    lines.append("%s %s %s %s\t%s\n" % (atom,res,chain,res_name,str(siss)))

    lines = sorted(lines,key=lambda x:int(x.split()[0]))

    lines = ["Residue\tSphereCon value\n"] + lines

    f = open(output_file,'w')
    f.write("".join(lines))
    f.close()

def calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha,atom_mode=False,manu_thresh=None,manu_angle_thresh=None,double_unknown_mode=False,manu_parameters=None,unknown_parameters=None,dist_matrix_only=False,res_name_map=None,sparsity=None):
    global threshs
    global angle_threshs
    global threshs_alpha
    global angle_threshs_alpha
    global unknown_thresh
    global unknown_angle_thresh
    global dist_matrix_threshs
    siss_map = {}

    if sparsity == None:
        sparsity = 'None'

    if unknown_parameters != None:
        (unknown_thresh,unknown_angle_thresh) = unknown_parameters
    if manu_parameters != None:
        [threshs,angle_threshs] = manu_parameters
        [threshs_alpha,angle_threshs_alpha] = manu_parameters
        dist_matrix_threshs = {sparsity:(manu_parameters[0])}
    if c_alpha:

        if dist_matrix_only:
            coordinate_map = res_name_map
        for chain in target_residues:
            siss_map[chain] = {}
            for res in target_residues[chain]:
                if not dist_matrix_only:
                    res_name = coordinate_map[chain][res][0]
                else:
                    res_name = res_name_map[chain][res]

                if manu_thresh == None:                
                    if dist_matrix_only:
                        if res_name in dist_matrix_threshs[sparsity]:
                            thresh = dist_matrix_threshs[sparsity][res_name]
                        else:
                            thresh = dist_matrix_threshs[sparsity]['UNK']
                    else:
                        thresh = threshs_alpha[res_name]
                else:
                    thresh = manu_thresh

                if manu_angle_thresh == None:
                    angle_thresh = angle_threshs_alpha[res_name]
                else:
                    angle_thresh = manu_angle_thresh

                if thresh == None or double_unknown_mode:
                    thresh = unknown_thresh
                if angle_thresh == None or double_unknown_mode:
                    angle_thresh = unknown_angle_thresh

                if dist_matrix_only:
                    angle_thresh = -1.0

                if not dist_matrix_only:
                    (atomname,x,y,z) = list(coordinate_map[chain][res][1].values())[0]
                    c_alpha_1 = numpy.array([x,y,z])
                    centroid_1 = centroid_map[chain][res][0]

                siss = 0.0
                for chain2 in coordinate_map:
                    for res_2 in coordinate_map[chain2]:
                        if chain == chain2 and res == res_2:
                            continue
                        if not dist_matrix_only:
                            res_name_2 = coordinate_map[chain2][res_2][0]
                        else:
                            res_name_2 = res_name_map[chain2][res_2]
                        if (chain,res,chain2,res_2) in dist_matrix:
                            d = dist_matrix[(chain,res,chain2,res_2)]
                        else:
                            d = None
                            #raise NameError('Did not find residue pair in distance matrix for: %s and %s' % (res,res_2))
                        if d == None: #This can happen for a sparse distance matrix
                            continue
                        rad = radii_map[res_name_2]
                        if double_unknown_mode:
                            rad = unknown_rad
                        if d - rad <= thresh:
                            if angle_thresh > -1.0: 
                                centroid_2 = centroid_map[chain2][res_2][0]
                                
                                angle = getCosAngle(centroid_1-c_alpha_1,centroid_2-c_alpha_1)
                                if angle == None:    
                                    siss += sphere_intersection(thresh,rad,d)
                                elif angle_thresh <= angle:
                                    siss += sphere_intersection(thresh,rad,d)
                            else:
                                siss += sphere_intersection(thresh,rad,d)

                rad = radii_map[res_name]
                if double_unknown_mode:
                    rad = unknown_rad
                siss = sissCorrectionVol(siss,thresh,angle_thresh,rad)
                siss_map[chain][res] = siss
    elif not atom_mode:
        for chain in target_residues:
            siss_map[chain] = {}
            for res in target_residues[chain]:
                if not res in coordinate_map[chain]:
                    continue
                res_name = coordinate_map[chain][res][0]
                atomlist = coordinate_map[chain][res][1]
                if manu_thresh == None:
                    thresh = threshs[res_name]
                else:
                    thresh = manu_thresh
                if manu_angle_thresh == None:
                    angle_thresh = angle_threshs[res_name]
                else:
                    angle_thresh = manu_angle_thresh

                if thresh == None or double_unknown_mode:
                    thresh = unknown_thresh
                if angle_thresh == None or double_unknown_mode:
                    angle_thresh = unknown_angle_thresh
                centroid_1 = centroid_map[chain][res]

                gly_c = [None]
                gly_n = [None]
                c_alpha_1 = [None]
                c_beta = [None]
                for atom in atomlist:
                    (atomname,x,y,z) = atomlist[atom]
                    if atomname == 'CA':
                        c_alpha_1 = numpy.array([x,y,z])
                    if atomname == 'CB':
                        c_beta = numpy.array([x,y,z])
                    if atomname == 'N':
                        gly_n = numpy.array([x,y,z])
                    if atomname == 'C':
                        gly_c = numpy.array([x,y,z])

                siss = 0.0
                error_flag = False
                for chain2,res_2 in dist_matrix[(chain,res)]:
                    for atom in dist_matrix[(chain,res)][(chain2,res_2)]:
                        
                        (d,atomname,x,y,z) = dist_matrix[(chain,res)][(chain2,res_2)][atom]
                        if d - vdw_radius[atomname[0]] <= thresh:
                            if angle_thresh > -1.0:
                                coord = numpy.array([x,y,z])

                                if res_name != 'GLY':
                                    if c_beta[0] == None or c_alpha_1[0] == None:
                                        angle = None
                                    else:
                                        angle = getCosAngle(c_beta-c_alpha_1,coord-c_alpha_1)

                                    if angle == None:
                                        error_flag = True
                                        #This happens for resiudes (beside Glycin), where only the C-Alpha atom is given. (Note: This is the full atom case)
                                        #The solution is to handle it as Glycin
                                        
                                        if gly_c[0] == None or gly_n[0] == None or c_alpha_1[0] == None:
                                            angle = 1.0
                                        else:
                                            gly_vec = gly_vector(gly_n,gly_c,c_alpha_1)
                                            angle = getCosAngle(gly_vec,coord-centroid_1)
                                    
                                else:
                                    if gly_c[0] == None or gly_n[0] == None or c_alpha_1[0] == None:
                                        error_flag = True
                                        angle = 1.0
                                    else:    
                                        gly_vec = gly_vector(gly_n,gly_c,c_alpha_1)
                                        angle = getCosAngle(gly_vec,coord-centroid_1)
                                
                                if angle == None or angle_thresh <= angle:
                                    siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)

                            else:
                                #If the full sphere is taken, there is no need for calculating the angle
                                siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)

                        #print(res_2,atomname,d,angle,siss)
                #print('\n')
                siss = sissCorrectionVol(siss,thresh,angle_thresh,radii_map[res_name])
                #print(siss,'\n')
                siss_map[chain][res] = siss

    else:
        for chain in target_residues:
            siss_map[chain] = {}
            for res in target_residues[chain]:
                siss_map[chain][res] = {}
                if not res in coordinate_map[chain]:
                    continue
                res_name = coordinate_map[chain][res][0]
                if manu_thresh == None:
                    thresh = threshs[res_name]
                else:
                    thresh = manu_thresh

                for chain2,res_2 in dist_matrix[(chain,res)]:
                    for atom1 in dist_matrix[(chain,res)][(chain2,res_2)]:
                        siss = 0.0
                        for atom2 in dist_matrix[(chain,res)][(chain2,res_2)][atom1]:
                            (d,atomname2,x2,y2,z2) = dist_matrix[(chain,res)][(chain2,res_2)][atom1][atom2]
                            
                            if d - vdw_radius[atomname2[0]] <= thresh:
                                #If the full sphere is taken, there is no need for calculating the angle
                                siss += sphere_intersection(thresh,vdw_radius[atomname2[0]],d)

                        siss = sissCorrectionVol(siss,thresh,-1.0)
                        siss_map[chain][res][atom1] = siss

    return siss_map

def calcAASiss(coordinate_map,target_residues,manu_thresh=None):
    siss_map = {}    
    
    for res in target_residues:
        siss = 0.0
        res_name = coordinate_map[res][0]
        atomlist = coordinate_map[res][1]
        
        if manu_thresh == None:
            thresh = threshs[res_name]
        else:
            thresh = manu_thresh
        
        for atom in atomlist:
            atom_siss = 0.0
            (atomname,x,y,z) = atomlist[atom]    
            for res_2 in coordinate_map:
                if res != res_2:
                    atomlist_2 = coordinate_map[res_2][1]
                    for atom_2 in atomlist_2:
                        (atomname_2,x_2,y_2,z_2) = atomlist_2[atom_2]
                        diff = numpy.array([x_2,y_2,z_2]) - numpy.array([x,y,z])
                        d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                        #If one atom is further than 10A+threshhold, then ignore the whole residue
                        if d - 15.0 > thresh:
                            break
                        if d - vdw_radius[atomname_2[0]] <= thresh:                    
                            atom_siss += sphere_intersection(thresh,vdw_radius[atomname_2[0]],d)

            siss += sissCorrectionVol(atom_siss,thresh,-1.0)
        siss_map[res] = siss/float(len(atomlist))       

    return siss_map

def siss(input_file=None,output_file=None,target_residues=None,chains=None,atom_mode=False,c_alpha=False,double_unknown_mode=False,dist_matrix_only=False,dist_matrix=None,res_name_map=None):
    
    cwd = os.getcwd()
    if output_file == None:
        if input_file == None:
            output_file = "%s/spherecon_output.tsv" % cwd
        else:
            output_file = "%s_spherecon.tsv" % input_file

    if not dist_matrix_only:
        if input_file == None:
            raise NameError("Error: no input file given.")
        else:
            coordinate_map,protein_centroid = parsePDB(input_file,chains,c_alpha)

        if target_residues == None:
            target_residues = {}
            for chain in coordinate_map:
                target_residues[chain] = list(coordinate_map[chain].keys())

        centroid_map = calcCentroidMap(coordinate_map,target_residues,c_alpha,double_unknown_mode=double_unknown_mode)

        dist_matrix = calcDistMatrix(coordinate_map,centroid_map,target_residues,c_alpha,atom_mode=atom_mode)
        sparsity = None
    else:
        if dist_matrix == None:
            dist_matrix,res_name_map,seq,sparsity = parseDM(input_file)
            #raise NameError("Error: Distance Matrix mode without given distance matrix.")
        if res_name_map == None:
            raise NameError("Error: Distance Matrix mode without given Residue Name Map.")
        if not c_alpha:
            print('Warning: Distance Matrix mode only possible in C-alpha only mode => C-alpha mode automatically activated.')
            c_alpha = True
        coordinate_map = None
        centroid_map = None

        if target_residues == None:
            target_residues = {'A':list(res_name_map.keys())}

    siss_map = calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha,atom_mode=atom_mode,double_unknown_mode=double_unknown_mode,dist_matrix_only=dist_matrix_only,res_name_map=res_name_map,sparsity=sparsity)
    
    produceOutput(siss_map,coordinate_map,output_file,res_name_map,atom_mode=atom_mode)
    
if __name__ == "__main__":
    
    argv = sys.argv[1:]
    helptext = """
\tSphereCon

This tool calculates a measure for relative solvent accessible area for single residues or whole amino acid chains.
Usage:
spherecon.py -i /Path/To/Input/File [-o /Path/To/Output/File] [-c chain] [-r residues] [--ca] [--bb] [--dm]

-i:\tPath to an input file in PDB (Protein Data Bank) file format or distance matrix file in case of --dm.

-o:\tPath to the output file produced as tab separated text file.
\tDefault: *InputFile*_spherecon.tsv

-c:\tChain identifier of the amino acid chain, or a list of Chain identifier separated by a ',', which should be analysed denoted as the chain identifiers of the ATOM records in the PDB file.
\tDefault: The first chain found in the input file
\tExamples: A | A,B | C,B,z,2

-r:\tList of residue identifiers for all residues for which the SphereCon value should be computed. If not given, the SphereCon values for all residues are computed.
\tThe residue identifiers denote as the residues identifiers of the ATOM records in the PDB file.
\tExamples: 234,78,368 | 17 | 34,35,36,37

--ca:\tC alpha version of SphereCon. Needs only coordinates of the C alpha atoms and their amino acid type.

--bb:\tBackbone only version of SphereCon. Needs only coordinates of the C alpha atoms.

--dm:\tDistance matrix version of SphereCon. Takes a distance matrix as input instead of a PDB file.
\tFile for distance matrix has to be in the following format:
\t\tSEQ [amino acid sequence]
\t\t[residue nr 1] [residue nr 2] [distance value in angstrom]
\t\t...
\tExample:
\t\tSEQ TGLRWGGDGIVQIVANNAIVGGWNSTDIFTEAGKHITSN
\t\t4 9 12.641
\t\t5 8 10.656
\t\t5 9 12.616
\t\t5 10 13.342
\t\t5 37 11.964
"""
    if len(argv) == 0:
        print(helptext)
        sys.exit()

    input_file = None
    output_file = None
    target_residues = None
    chains = None
    c_alpha = False
    backbone = False
    atom_mode = False

    dist_matrix_only = False

    try:
        opts,args = getopt.getopt(argv,"ahr:o:i:c:",['ca','bb','dm'])
    except getopt.GetoptError:
        print("spherecon.py -h")
        sys.exit(2)
    for opt,arg in opts:
        if opt == '-h':
            print(helptext)
            sys.exit()
        elif opt == '-a':
            atom_mode = True
        elif opt == "-i":
            input_file = arg
        elif opt == '-o':
            output_file = arg
        elif opt == '-r':
            target_residues = arg.split(',')
        elif opt == '-c':
            chain_str = arg
            chains = set(chain_str.split(','))
        elif opt == '--ca':
            c_alpha = True
        elif opt == '--bb':
            backbone = True
            c_alpha = True
        elif opt == '--dm':
            dist_matrix_only = True
            c_alpha = True

    siss(input_file=input_file,output_file=output_file,target_residues=target_residues,chains=chains,atom_mode=atom_mode,c_alpha=c_alpha,double_unknown_mode=backbone,dist_matrix_only=dist_matrix_only,dist_matrix=None,res_name_map=None)

