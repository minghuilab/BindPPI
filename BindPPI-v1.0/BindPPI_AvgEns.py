#!/usr/bin/env python
# coding: utf-8

import random
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor
from torch.utils.data import DataLoader 
from models import *
import sys
esm_path = './esm/'
sys.path.append(esm_path)
import esm
import json
import argparse
import os, re
from collections import defaultdict,Counter,deque
from string import ascii_uppercase
from string import ascii_lowercase
import subprocess, sys, getopt
import joblib
import linecache
import warnings 
warnings.filterwarnings('ignore')


ascii_cases = ascii_uppercase + ascii_lowercase




map_3_to_1 = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
              "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
              "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
              "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
              "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
              "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
              "PTR": "X", "XLE": "X", "XAA": "X", "HSD": "H"}
# standard state amino acid suface areas.REF: Hydrophobicity of amino acid residues in globular proteins
map_surface = {'A':118.1,'R':256.0,'N':165.5,'D':158.7,'C':146.1,'Q':193.2,'E':186.2,'G':88.1,'H':202.5,'I':181.0,'L':193.1,'K':225.8,'M':203.4,'F':222.8,'P':146.8,'S':129.8,'T':152.5,'W':266.3,'Y':236.8,'V':164.5}
aa_cpa = {'K':'C', 'R':'C', 'D':'C', 'E':'C',###charge
          'H':'P', 'S':'P', 'T':'P', 'N':'P', 'Q':'P', 'C':'P', 'Y':'P',###polar
          'A':'A', 'G':'A', 'I':'A', 'L':'A', 'F':'A', 'P':'A', 'V':'A', 'M':'A', 'W':'A'}###apolar

jobid = ''
myopts, args = getopt.getopt(sys.argv[1:], "i:")
for o, a in myopts:
    if o == '-i':
        jobid = a
    else:
        print ("Usage: %s -i jobid" % sys.argv[0])

workdir = Your working directory # '/data/jiang/tools/bindppi_test/' ##need change, not contain capital letters
pathinput = path for inputfiles directory # '/data/jiang/tools/bindppi_test/inputfiles' ##need change, not contain capital letters
pathvmd = path for VMD software # '/usr/local/bin/vmd' ##need change
pathcharmm = path for CHARMM software # '/usr/local/bin/charmm' ##need change
pathnamd2 = path for NAMD software # '/usr/local/bin/namd2' ##need change
pathdssp = path for DSSP software # '/usr/local/bin/mkdssp' ##need change
pathipot = path for iPot software # '/data/jiang/tools/iPot/' ##need change
pathmcvol = path for McVol software #'/data/jiang/tools/McVol.rev/' ##need change
pathpdb2pqr30 = path for PDB2PQR software #'/data/jiang/anaconda3/bin/pdb2pqr30' ##need change
pathprovean = path for PROVEAN software # '/usr/local/bin/provean.sh' ##need change


path_model = pathinput+'/RF13.pkl'
jobpath = workdir + jobid
pathoutput = workdir + jobid + '_out'
os.system("mkdir %s" % pathoutput)
in_file = jobpath + '/' + jobid + '.input'


need_feas = ['L_p1/SA_p1','L_p2/SA_p2','P_Helix','P_Sheet','SA_p1','SA_p2','DE_elec','P_CS','E_ce','N_AP','S_AA','P_Charge','N_Water']

H_flag = 'withoutH'
distance_cutoff = '6'
suffix = H_flag+'_'+distance_cutoff


# # structure minimization

def ProPDB1():
    pdball = []
    with open(in_file, 'r') as f:
        f.readline()
        for line in f:
            ff = line.split('\t')
            pdb = ff[0]  # 1AUD.pdb or 1AUD.pdb1
            if pdb not in pdball:
                pdball.append(pdb)
                if pdb.split(".")[1] == "pdb1":
                    ST_play = False
                    ffpdb = open(pathoutput + '/' + pdb.split(".")[0].lower() + '_p.pdb', 'w')
                    fpdb = open(jobpath + '/' + pdb, 'r')
                    for linepdb in fpdb:
                        if linepdb[0:5] == "MODEL":
                            CountModel = linepdb.split()[1]
                            ST_play = True
                            continue
                        if ST_play:
                            if linepdb[:4] == "ATOM":
                                ffpdb.write("%s                 %s\n" % (
                                linepdb[0:55].strip('\r\n'), str(linepdb[21:22]) + '_' + str(CountModel)))
                else:
                    ffpdb = open(pathoutput + '/' + pdb.split(".")[0].lower() + '_p.pdb', 'w')  # 1aud_p.pdb
                    fpdb = open(jobpath + '/' + pdb, 'r')
                    for linepdb in fpdb:
                        line_list = re.split(r'\s+', linepdb)
                        if (line_list[0] == 'MODEL') and (line_list[1] == '2'):
                            break
                        if linepdb[:4] == 'ATOM':
                            ffpdb.write("%s                 %s\n" % (linepdb[0:55].strip('\r\n'), str(linepdb[21:22]) + '_' + str(1)))
            else:
                continue
    ffpdb.close()
    fpdb.close()

def del_unknown_incomplete():
    residues = [k for k, v in map_3_to_1.items() if v!= 'X']
    f = open(in_file, 'r')
    f.readline()
    for line in f:
        ff = line.split("\t")
        pdbid = ff[0].split('.')[0].lower()
        pdbfile = open('{}/{}_p.pdb'.format(pathoutput,pdbid)).readlines()
        # delete unknow resdues
        pdbfile_del_unknown = [i for i in pdbfile if i[17:20] in residues]
        # delete incomplete residues
        final_row = pdbfile_del_unknown[-1]
        last = ''
        above = []
        allresidues = []
        for row in pdbfile_del_unknown:
            if row[17:26] == last and row == final_row:            # when read final row，append it if equal to last
                if row [16] in [' ','A']:
                    above.append(row)
                atoms = [i[13:16].strip() for i in above if (i[16]==' ')|(i[16]=='A')] ##2j12:55.633  23.714 -12.966
                if set(['C','N','O','CA']).issubset(set(atoms)):
                    allresidues.append(above)
            elif row[17:26] == last and row != final_row:            # when read same residue, but not last row
                if row [16] in [' ','A']:
                    above.append(row)
            else:                            # when read different residue
                if len(above)>=4:
                    atoms = [i[13:16].strip() for i in above if (i[16]==' ')|(i[16]=='A')] ##2j12:55.633  23.714 -12.966
                    if set(['C','N','O','CA']).issubset(set(atoms)):
                        allresidues.append(above)
                above = [row]
            last = row[17:26]
        # write out
        with open('{}/{}_p_test.pdb'.format(pathoutput,pdbid),'w') as fw:
            fw.write(''.join([j for i in allresidues for j in i]))
        break
    os.system('mv {}/{}_p_test.pdb {}/{}_p.pdb'.format(pathoutput,pdbid, pathoutput,pdbid))
    f.close()

def splitchain():
    pdball = []
    f = open(in_file, 'r')
    f.readline()
    for line in f:
        ff = line.split("\t")
        pdbfile = ff[0]
        pdb = pdbfile.split(".")[0].upper()
        partner1 = ff[1].split(".")
        partner2 = ff[2].strip().split(".")
        if pdbfile not in pdball:
            for chains in list(partner1 + partner2):
                f1 = open(pathoutput + '/' + pdb.lower() + '_p.pdb', 'r')
                fw = open(pathoutput + '/' + pdb + '_' + chains + '.pdb', 'w')
                for line1 in f1:
                    if chains == line1[72:].strip(): 
                        fw.write(line1)
                f1.close()
                fw.close()
            pdball.append(pdbfile)
        else:
            continue
    f.close()

def CleanPdb():
    first_line = open(in_file, 'r').readlines()[0][:-1]
    fw = open(in_file + ".cleaned", "w")
    fw.write("%s\t%s\t%s\t%s\n" % (first_line, "PDBid", "NewPartner1", "NewPartner2"))

    second_line = open(in_file, 'r').readlines()[1][:-1]
    ff = second_line.split("\t")
    pdb = ff[0].split(".")[0].upper()
    partner1 = ff[1].split(".")
    partner2 = ff[2].strip().split(".")

    mapchainarray = []
    counti = 0
    for chains in list(partner1 + partner2):
        cc = (chains, ascii_cases[counti])
        mapchainarray.append(cc)
        counti += 1
        mapchaindict = dict(iter(mapchainarray))

    newpartner1 = ''
    for chains in list(partner1):
        newpartner1 += mapchaindict[chains]
    newpartner2 = ''
    for chains in list(partner2):
        newpartner2 += mapchaindict[chains]
    
    fw.write("%s\t%s\t%s\t%s\n" % (second_line, pdb, newpartner1, newpartner2))
    
    countchain = 1
    for chains in list(partner1 + partner2):
        fwpdb = open(pathoutput + "/" + pdb + "_CH" + str(countchain) + ".pdb", "w")
        countchain += 1
        count = 1
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        resname = fpdb.readlines()[0][17:20].strip()
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        resnum = fpdb.readlines()[0][22:27].strip()
        line1 = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r").readlines()[0]
        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        atomname = fpdb.readlines()[0][13:16].strip()
        fwpdb.write("%s %s%s   %s %s" % (line1[0:16], line1[17:21], mapchaindict[chains], str(count), line1[27:]))

        fpdb = open(pathoutput + "/" + pdb + "_" + chains + ".pdb", "r")
        for linepdb in fpdb.readlines()[1:]:
            resnamepdb = linepdb[17:20].strip()
            resnumpdb = linepdb[22:27].strip()
            atomnamepdb = linepdb[13:16].strip()
            if resnamepdb == resname and resnumpdb == resnum:
                if atomnamepdb != atomname:
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 10 and count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 100 and count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 1000 and count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                else:
                    continue
            else:
                if atomnamepdb != atomname:
                    count += 1
                    if count < 10:
                        fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 10 and count < 100:
                        fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 100 and count < 1000:
                        fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    if count >= 1000 and count < 10000:
                        fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
#                     print("DEBUG: how often we go in this branch", resnamepdb, resname, resnumpdb, resnum, atomnamepdb, atomname)
                else:
                    continue

            resname = linepdb[17:20].strip()
            resnum = linepdb[22:27].strip()
            atomname = linepdb[13:16].strip()

        fpdb.close()
        fwpdb.close()
    fw.close()

def wtpdb():
    second_line = linecache.getline(in_file+'.cleaned',2)[:-1]
    p1 = list(second_line.split('\t')[4])
    p2 = list(second_line.split('\t')[5])
    pdb = second_line.split('\t')[3]
    pdball = []
    if pdb not in pdball:
        subprocess.getoutput('cat %s/%s_CH*.pdb > %s/%s.pdb' % (pathoutput, pdb, pathoutput, pdb))
        pdball.append(pdb)
    return second_line,p1,p2,pdb

def vmd_wt():
    pdball = []
    template = open(pathinput + '/vmd.pgn').read()
    NumChain = int(len(p1 + p2))
    if pdb not in pdball:
        vmd_pdb = template.replace('protname', pdb).replace('NumChain', str(NumChain)).replace('pathinput', pathinput).replace('pathoutput', pathoutput)
        with open(pathoutput + '/vmd_' + pdb + '.pgn', 'w') as fw:
            fw.write(vmd_pdb)
        subprocess.getoutput('%s -dispdev text -e %s/vmd_%s.pgn' % (pathvmd, pathoutput, pdb))
        pdball.append(pdb)

def charmmfile_wt():
    pdball = []
    NumP = int(len(p1 + p2))
    countp = 1
    pdb_low = pdb.lower()
    if pdb_low not in pdball:
        countp = 1
        for chains in (list(p1) + list(p2)):
            os.system('grep "^.\{21\}%s" %s/%s_vmd.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb.upper(), pathoutput, pdb_low, 'ch' + str(countp)))
            ffpdb = open(pathoutput + '/' + pdb_low + '_ch' + str(countp) + '.pdb', 'a')
            ffpdb.write("%s" % ('END'))
            ffpdb.close()
            os.system('%s <%s/setup.inp path=%s protname=%s I=%s pathpara=%s> %s/setup_%s.out' % (pathcharmm,  pathinput, pathoutput, pdb_low, countp, pathinput, pathoutput, pdb_low + '_ch' + str(countp)))
            countp += 1
        os.system('%s <%s/append.inp path=%s protname=%s NumP=%s pathpara=%s> %s/append_%s.out' % (pathcharmm, pathinput, pathoutput, pdb_low, NumP, pathinput, pathoutput, pdb_low))
        os.system('%s <%s/charmm_to_namd.inp path=%s protname=%s pathpara=%s> %s/charmm_to_namd_%s.out' % (pathcharmm, pathinput, pathoutput, pdb_low, pathinput, pathoutput, pdb_low))
        pdball.append(pdb_low)

def runrest():
    pdball = []
    template = open(pathinput + '/rest.pgn').read()
    pdb_low = pdb.lower()
    if pdb_low not in pdball:
        rest_wild = template.replace('protname', pdb_low).replace('pathoutput', pathoutput)
        with open(pathoutput + '/rest_' + pdb_low + '.pgn', 'w') as fw:
            fw.write(rest_wild)
        subprocess.getoutput('%s -dispdev text -e %s/rest_%s.pgn' % (pathvmd, pathoutput, pdb_low))
        pdball.append(pdb_low)

def minimization():
    pdball = []
    runstep = "100"
    pdb_low = pdb.lower()
    template = open(pathinput + '/min.conf').read()
    if pdb_low not in pdball:
        min_wild = template.replace('protname', pdb_low).replace('pathoutput', pathoutput).replace('runstep', runstep).replace('pathinput', pathinput)
        with open(pathoutput + '/min_' + pdb_low + '.conf', 'w') as fw:
            fw.write(min_wild)
        subprocess.getoutput('%s %s/min_%s.conf' % (pathnamd2, pathoutput, pdb_low))
        pdball.append(pdb_low)

def frames():
    pdball = []
    pdb_low = pdb.lower()
    startframe = 100
    if pdb_low not in pdball:
        os.system('%s <%s/frame.inp path=%s protname=%s startframe=%s pathinput=%s> %s/frame_%s.out' % (pathcharmm, pathinput, pathoutput, pdb_low, startframe, pathinput, pathoutput, pdb_low))
        pdball.append(pdb_low)

def min_add_chain():
    pdb_low = pdb.lower()
    f = open(pathoutput+'/'+pdb_low+'_min.pdb','r')
    ff = open(pathoutput+'/'+pdb+'_min.pdb','w')
    chain_dict = dict(zip(['CH'+str(i) for i in range(1,len(p1+p2)+1)],p1+p2))
    for line in f.readlines():
        if line[:4]=='ATOM':
            suff_chain = line[70:-1].strip()
            new_line = line[:21]+chain_dict[suff_chain]+line[22:]
            need_line = new_line.replace('HSD','HIS').replace('OT1','O  ').replace('OT2','OXT')
            ff.write(need_line)
    f.close()
    ff.close()


# # features

# ## tran_pdb_to_dict

def tran_pdb_to_dict():
    pdb_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(np.array)))
    f = open(pathoutput+'/'+pdb+'_min.pdb','r')
    for line in f.readlines():
        if line[:4]=='ATOM':
            atom_number = line[4:11].strip()
            atom = line[11:17].strip()
            residue = map_3_to_1[line[17:20].strip()]
            chain = line[21]
            residue_number = line[22:27].strip()
            res = residue+chain+residue_number
            x_coord = float(line[27:38].strip())
            y_coord = float(line[38:46].strip())
            z_coord = float(line[46:54].strip())
            if line[11:16].strip()[0]!='H':
                pdb_dict[chain][res][atom_number+'_'+atom] = np.array([x_coord,y_coord,z_coord])
    f.close()
    return pdb_dict


# ## energy 

def energyfile():
    """
    produce energy inputfile
        """
    partnerall = []
    f = open(in_file + ".cleaned", 'r')
    f.readline()
    template = open(pathinput + '/energy.inp').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[3].lower()
        partner1 = ff[4]
        partner2 = ff[5]
        NumP1 = int(len(partner1))
        NumP2 = int(len(partner2))
        p1 = 'ch1'
        p2 = 'ch' + str(NumP1 + 1)
        for countp1 in range(2, NumP1 + 1):
            p1 += ' -\n .or. segid ch' + str(countp1)
        for countp2 in range(NumP1 + 2, NumP1 + NumP2 + 1):
            p2 += ' -\n .or. segid ch' + str(countp2)
        mp = p1
        op = p2
        partner = 'p1'
        partnerall = []
        if partner not in partnerall:
            energy_wild = template.replace('p1', mp).replace('p2', op).replace('pathinput', pathinput)
            fwild = open(pathoutput + "/energy_" + pdb+'_'+partner +".inp", "w")
            fwild.write(energy_wild)
            partnerall.append(partner)
            fwild.close()
    f.close()

def energy():
    """
        energy calculation
        """
    partnerall = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.readline()
    for line in f:
        ff = line.split("\t")
        pdb = ff[3].lower()
        partner1 = ff[4]
        partner2 = ff[5]
        partner = 'p1'
        if partner not in partnerall:
            os.system('%s <%s/energy_%s.inp path=%s protname=%s > %s/energy_%s.out' % (pathcharmm, pathoutput, pdb+'_'+partner, pathoutput, pdb, pathoutput, pdb+'_'+partner))
            partnerall.append(partner)
    f.close()


def get_energy_asa_len():
    elec_wild = float(os.popen('grep \'ELECENER <\' %s/energy_%s.out | cut -c 25-100' % (pathoutput, pdb.lower()+'_p1')).read().strip()[1:-1])
    msmp_wild = float(os.popen('grep \'MSA <\' %s/energy_%s.out | cut -c 20-100' % (pathoutput, pdb.lower()+'_p1')).read().strip()[1:-1])
    msop_wild = float(os.popen('grep \'MSC <\' %s/energy_%s.out | cut -c 20-100' % (pathoutput, pdb.lower()+'_p1')).read().strip()[1:-1])
    len_p1 = sum([len(pdb_dict[chain]) for chain in p1])
    len_p2 = sum([len(pdb_dict[chain]) for chain in p2])
    lenp1_ASAP1 = len_p1/msmp_wild if len_p1>=len_p2 else len_p2/msop_wild
    lenp2_ASAP2 = len_p2/msop_wild if len_p1>=len_p2 else len_p1/msmp_wild
    msp1wt = msmp_wild if len_p1>=len_p2 else msop_wild
    msp2wt = msop_wild if len_p1>=len_p2 else msmp_wild
    with open(pathoutput+'/'+pdb+'_energy_asa_len.'+suffix,'w') as ff:
        ff.write('DE_elec\tL_p1/SA_p1\tL_p2/SA_p2\tSA_p1\tSA_p2\n')
        ff.write('%s\t%s\t%s\t%s\t%s\n'%(elec_wild,lenp1_ASAP1,lenp2_ASAP2,msp1wt,msp2wt))


# ## distance

def distance_contact_resid():
    with open(in_file + ".cleaned", 'r') as f:
        f.readline()
        for line in f:
            ff = line.split("\t")
            pdb = ff[3].lower()  # 1aay
            partner1 = ff[4]  # cleaned
            partner2 = ff[5]
            NumP1 = int(len(partner1))
            NumP2 = int(len(partner2))
            p1 = ['ch{}'.format(i) for i in range(1, NumP1 + 1)]
            p1 = ' -\n .or. segid '.join(p1)
            p2 = ['ch{}'.format(i) for i in range(NumP1 + 1, NumP1 + NumP2 + 1)]
            p2 = ' -\n .or. segid '.join(p2)
            mp = p1
            op = p2
            for cutoff in [distance_cutoff]:
                cuthb = cutoff.split('.')[0]
                # wild type
                with open('{}/pdb_distance_withH.inp'.format(pathinput), 'r') as f1:
                    template = f1.read()
                result = template.replace('p1', mp).replace('p2', op).replace('cutoff', cutoff)
                with open('{}/distance_{}_{}.inp'.format(pathoutput, pdb, cuthb), 'w') as fwild:
                    fwild.write(result)
                os.system('%s <%s/distance_%s_%s.inp path=%s protname=%s pathpara=%s > %s/distance_%s_%s.out' % (pathcharmm, pathoutput, pdb, cuthb, pathoutput, pdb, pathinput, pathoutput, pdb, cuthb))

def distance_contact_extract():
    chain_dict = dict(zip(['CH'+str(i) for i in range(1,len(p1+p2)+1)],p1+p2))
    distance_count = defaultdict(lambda: defaultdict())
    distance_count[H_flag][distance_cutoff]=[]
    cutoff = float(distance_cutoff)
    with open(pathoutput + '/' + 'distance_' + pdb.lower() + '_'+distance_cutoff+'.out', 'r') as f:
        for line in f.readlines():
            if re.match(r'\s+\d+ CH', line) or re.match(r'\d+ CH', line):
                res1 = map_3_to_1[line[11:14]]+chain_dict[line[6:9]]+line[14:21].strip()
                atom1 = line[21:25].strip()
                res2 = map_3_to_1[line[36:39]]+chain_dict[line[31:34]]+line[39:46].strip()
                atom2 = line[46:50].strip()
                res1,atom1,res2,atom2
                distance = float(line[53:62].strip())
                if distance < cutoff:
                    contact = (res1,res2)
                    if H_flag == 'withoutH':
                        if (contact not in distance_count[H_flag][distance_cutoff])&(atom1[0]!='H')&(atom2[0]!='H'):
                            distance_count[H_flag][distance_cutoff].append(contact)
                    else:
                        if (contact not in distance_count[H_flag][distance_cutoff]):
                            distance_count[H_flag][distance_cutoff].append(contact)
    return distance_count


# ## contact pairs

def get_contacts():
    cutoff = int(suffix.split('_')[1])
    chain_dict = dict(zip(['CH'+str(i) for i in range(1,len(p1+p2)+1)],p1+p2))
    with open(pathoutput + '/' + 'distance_' + pdb.lower() + '_'+distance_cutoff+'.out', 'r') as f,open(pathoutput+'/'+pdb+'_contacts.'+suffix,'w') as ff:
        atom_pairs = []
        for line in f.readlines():
            if re.match(r'\s+\d+ CH', line) or re.match(r'\d+ CH', line):
                res1 = map_3_to_1[line[11:14]]+chain_dict[line[6:9]]+line[14:21].strip()
                atom1 = line[21:25].strip()
                res2 = map_3_to_1[line[36:39]]+chain_dict[line[31:34]]+line[39:46].strip()
                atom2 = line[46:50].strip()
                distance = float(line[53:62].strip())
                if distance < cutoff:
                    contact_atom_pair = (res1,atom1,res2,atom2)
                    if (contact_atom_pair not in atom_pairs)&(atom1[0]!='H')&(atom2[0]!='H'):
                        atom_pairs.append(contact_atom_pair)
        ff.write('N_AP\n')
        ff.write('%s\n'%(len(atom_pairs)))


# ## dssp

def get_partner_pdb():
    f = open(pathoutput+'/'+pdb+'_min.pdb','r')
    f_p1 = open(pathoutput+'/'+pdb+'_min_p1.pdb','w')
    f_p2 = open(pathoutput+'/'+pdb+'_min_p2.pdb','w')
    for line in f.readlines():
        if line[:4]=='ATOM':
            if line[21] in p1:
                f_p1.write(line)
            if line[21] in p2:
                f_p2.write(line)
    f.close()
    f_p1.close()
    f_p2.close()

def run_dssp():
    pdball = []
    if pdb not in pdball:
        os.system('%s -i %s/%s_min.pdb -o %s/%s_min.dssp' % (pathdssp, pathoutput, pdb, pathoutput, pdb))
        os.system('%s -i %s/%s_min_p1.pdb -o %s/%s_min_p1.dssp' % (pathdssp, pathoutput, pdb, pathoutput, pdb))
        os.system('%s -i %s/%s_min_p2.pdb -o %s/%s_min_p2.dssp' % (pathdssp, pathoutput, pdb, pathoutput, pdb))
        pdball.append(pdb)
    pdb_dssp = defaultdict(str)
    f = open(pathoutput+'/'+pdb+'_min.dssp','r')
    ST_play = False
    for line in f.readlines():
        if line[:25]=='  #  RESIDUE AA STRUCTURE':
            ST_play = True
            continue
        if ST_play:
            if '!' not in line:
                res_dssp = line[13]+line[10:12].strip()+line[5:10].strip()
                res_sec = line[16]
                pdb_dssp[res_dssp] = res_sec
    f.close()
    return pdb_dssp

def generate_dssp_asa():
    pdb_dssp_asa = defaultdict(lambda:defaultdict(int))
    state_list = ['com','p1','p2']
    file_list = [pathoutput+'/'+pdb+'_min.dssp',pathoutput+'/'+pdb+'_min_p1.dssp',pathoutput+'/'+pdb+'_min_p2.dssp']
    state_file = dict(zip(state_list,file_list))
    for state in state_file:
        file = state_file[state]
        f = open(file,'r')
        ST_play = False
        for line in f.readlines():
            if line[:25]=='  #  RESIDUE AA STRUCTURE':
                ST_play = True
                continue
            if ST_play:
                if '!' not in line:
                    res_dssp = line[13]+line[10:12].strip()+line[5:10].strip()
                    res_acc = int(line[34:38].strip())
                    pdb_dssp_asa[state][res_dssp] = res_acc
        f.close()
    return pdb_dssp_asa

def get_secondary_structure():
    dssp_res = defaultdict(str)
    for res in residues:
        if pdb_dssp[res] ==' ':
            dssp_res[res] = 'C'
        else:
            dssp_res[res] = pdb_dssp[res]
    with open(pathoutput+'/'+pdb+'_inter_sec_stru.'+suffix,'w') as ff:
        ff.write('P_Sheet\tP_Helix\n')
        proportion_BE = len([dssp_res[i] for i in residues if dssp_res[i] in ['B','E']])/len(residues)# sheet
        proportion_HGI = len([dssp_res[i] for i in residues if dssp_res[i] in ['H','G','I']])/len(residues)# helix
        ff.write('%s\t%s\n'%(proportion_BE,proportion_HGI))


def cal_nis():
    surface_residues = []
    for res in pdb_dssp_asa['com']:
        rasac = pdb_dssp_asa['com'][res]/map_surface[res[0]]
        if rasac > 0.2: ###surface
            surface_residues.append(res)
    surface_residues1 = [i for i in surface_residues if i not in residues]
    f = open(pathoutput+'/'+pdb+'_nis.'+suffix,'w')
    f.write('P_Charge\n')
    dict_cpa_count = dict(Counter([aa_cpa[i[0]] for i in surface_residues1]))
    nis_list = []
    for cate in ['C']:
        cate_nis = ''
        if cate in dict_cpa_count:
            cate_nis = dict_cpa_count[cate]/len(surface_residues1)
        else:
            cate_nis = 0
        nis_list.append(cate_nis)
    f.write('%s\n'%(nis_list[0]))
    f.close()


# ## iPot

def cal_ipot():
    matrix_list = ['AACE167']
    dmax_list = [10.0]
    kmin_list = [5]
    pdb_low = pdb.lower()
    ipot_names = [matrix.lower()+'_d'+str(int(dmax))+'_k'+str(kmin) for matrix in matrix_list for dmax in dmax_list for kmin in kmin_list]
    receptor_path = pathoutput+'/'+pdb+'_min_p1.pdb'
    ligand_path = pathoutput+'/'+pdb+'_min_p2.pdb'
    pdb_ipot_ene_list = []
    matrix,dmax,kmin = 'AACE167','10.0','5'
    table = pathipot+'ipot_data/'+matrix+'/table.'+str(dmax)+'A_k'+str(kmin)
    ipot_ene = subprocess.getoutput('%s%s -r %s -l %s -t %s -d %s'%(
        pathipot,matrix.lower(),receptor_path,ligand_path,table,dmax)).split('=')[1].strip()
    with open(pathoutput+'/'+pdb+'_ipot.'+suffix,'w') as f:
        f.write('E_ce\n')
        f.write('%s\n'%(ipot_ene))


# ## McVol

# process mcvol outfile, and output the number, volume and waters for cavities and cleft.
def process_mcvol(mcvolfile):
    waters = {i.strip().split()[-1]:float(i.strip().split()[0]) for i in open(mcvolfile).readlines() if 'waters placed in' in i}  # cavity id : waters num
    cavities_cleft = [i for i in open(mcvolfile).readlines() if i.startswith('Checking ')]
    cavities = {i.strip().split()[2]:float(i.strip().split()[-1].replace('A^3','')) for i in cavities_cleft if i.startswith('Checking cavity ') and i.strip().split()[2] in waters.keys()}
    cleft = {i.strip().split()[2]:float(i.strip().split()[-1].replace('A^3','')) for i in cavities_cleft if i.startswith('Checking cleft ') and i.strip().split()[2] in waters.keys()}
    mcvol_cavity_water_whole = np.sum([waters[i] for i in cavities.keys()])
    return mcvol_cavity_water_whole

def run_McVol():
    f = open(pathoutput+'/'+pdb+'_min.pdb','r')
    ff = open(pathoutput+'/'+pdb+'_min_noH.pdb','w')
    for line in f.readlines():
        if line[11:17].strip()[0]!='H':
            ff.write(line)
    f.close()
    ff.close()
    os.system('%s --ff=CHARMM %s/%s_min_noH.pdb %s/%s_min.pqr'%(pathpdb2pqr30,pathoutput,pdb,pathoutput,pdb))
    if not os.path.exists(pathoutput+'/'+pdb+'_min.pqr'):
        print('pdb2pqr run error :%s'%(jobid))

    os.chdir(pathoutput)
    os.system('cp %sTest/all.setup %s/%s_min.setup'%(pathmcvol,pathoutput,pdb))
    cmd= '%ssrc/McVol %s_min > McVol_%s_whole.out'%(pathmcvol,pdb,pdb)
    cmdlist = deque([cmd])
    while cmdlist:
        cmd = cmdlist.popleft()
        process = subprocess.Popen(cmd,stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True, close_fds=True, start_new_session=True)
        try:
            output, unused_err = process.communicate(timeout=1200)
            print('mcvol run Success: '+jobid)
        except subprocess.TimeoutExpired:
            cmdlist.append(cmd)
            os.system('rm %s' % (pathoutput+'/McVol_'+pdb+'_whole.out'))
            try:
                os.killpg(process.pid, signal.SIGTERM)
                print('mcvol run Success Failed: '+jobid)
            except ProcessLookupError:
                process.kill()
                print('mcvol run Failed again: '+jobid)

        os.system('rm %s' % (pathoutput+'/surface_clusters_surf.pdb'))
        os.system('rm %s' % (pathoutput+'/surface_clusters_cav.pdb'))
        os.system('rm %s' % (pathoutput+'/test2.pdb'))
        os.system('rm -r %s' % (pathoutput+'/cav'))

    mcvol_out = pathoutput+'/McVol_'+pdb+'_whole.out'
    with open(pathoutput+'/'+pdb+'_McVol_whole.'+suffix,'w') as f1:
        f1.write('N_Water\n')
        mcvol_cavity_water_whole = process_mcvol(mcvol_out)
        f1.write('%s\n'%(mcvol_cavity_water_whole))


# ## provean

def gen_seq_var():
    for chain in pdb_dict:
        with open(pathoutput+'/'+pdb+'_'+chain+'.seq','w') as f_seq:
            chain_seq = ''.join([res[0] for res in pdb_dict[chain]])
            f_seq.write('> %s_%s\n'%(pdb,chain))
            f_seq.write(chain_seq)
    for chain in p1+p2:
        chain_residues = [i for i in residues if i[1]==chain]
        with open(pathoutput+'/'+pdb+'_'+chain+'.var','w') as fw:
            chain_seq = ''.join([res[0] for res in pdb_dict[chain]])
            for res in chain_residues:
                var = res[0]+res[2:]+'A'
                fw.write('%s\n'%(var))

def run_provean():
    chainall = []
    for chain in p1+p2:
        if chain not in chainall:
            chainall.append(chain)
            os.system('%s -q %s/%s_%s.seq -v %s/%s_%s.var > %s/provean_%s_%s.out --num_threads 30'%(
                pathprovean,pathoutput,pdb,chain,pathoutput,pdb,chain,pathoutput,pdb,chain))

def cal_proportion_cons_provean():
    provean_dict = {}
    for chain in p1+p2:
        provean_file = 'provean_'+pdb+'_'+chain+'.out'
        with open(pathoutput+'/'+provean_file,'r') as f,open(pathoutput+'/'+provean_file+'_gai','w') as ff:
            flag=False
            for line in f.readlines():
                if flag:
                    ff.write(line)
                if line[:11]=='# VARIATION':
                    flag = True
        provean_result = pd.read_csv(pathoutput+'/'+provean_file+'_gai',sep='\t',names=['mutation','provean'])
        provean_result['mutation1'] = provean_result['mutation'].str[0]+chain+provean_result['mutation'].str[1:]
        provean_dict1 = dict(zip(list(provean_result['mutation1']),list(provean_result['provean'])))
        provean_dict.update(provean_dict1)
    cons_residues = []
    for res in residues:
        if res+'A' in provean_dict:
            if provean_dict[res+'A'] <=-2.5:
                cons_residues.append(res)
    proportion_cons_provean = len(cons_residues)/len(residues)
    f = open(pathoutput+'/'+pdb+'_provean.'+suffix,'w')
    f.write('P_CS\n')
    f.write('%s\n'%(proportion_cons_provean))
    f.close()


# ## MITS020101
def cal_MITS020101():
    aaindex_1 = pd.read_csv(pathinput+'/aaindex_1',sep='\t')
    name = 'MITS020101'
    aaind1_dict = dict(zip(list(aaindex_1['aa']),list(aaindex_1[name])))
    name_sum = sum([aaind1_dict[res[0]] for res in residues])
    with open(pathoutput+'/'+pdb+'_aaindex1.'+suffix,'w') as f:
        f.write('S_AA\n')
        f.write('%s\n'%(name_sum))

#structure minimization
ProPDB1()
del_unknown_incomplete()
splitchain()
CleanPdb()
second_line,p1,p2,pdb = wtpdb()
vmd_wt()
charmmfile_wt()
runrest()
minimization()
frames()
min_add_chain()
#features
pdb_dict = tran_pdb_to_dict()
energyfile()
energy()
get_energy_asa_len()
distance_contact_resid()
distance_count = distance_contact_extract()
residues = sorted(set([j for i in distance_count[H_flag][distance_cutoff] for j in i]))
if len(residues)==0:
    print('No interfacial residue was observed between interaction partners (Interface residues are defined as those with inter-atomic distances less than 6Å between any heavy atoms of interacting partners).')
else:
    get_contacts()
    get_partner_pdb()
    pdb_dssp = run_dssp()
    pdb_dssp_asa = generate_dssp_asa()
    get_secondary_structure()
    cal_nis()
    cal_ipot()
    run_McVol()
    gen_seq_var()
    run_provean()
    cal_proportion_cons_provean()
    cal_MITS020101()


    embedding_file_path = pathoutput + '/'

    json_file_path = jobpath + '/' + jobid + '.json'

    with open(json_file_path, 'r') as fw:
        seq_dict = json.load(fw)

    model_name = "esm2_t36_3B_UR50D"
    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()#.to(device)
    def get_embed(datatmp):    
        batch_labels, batch_strs, batch_tokens = batch_converter(datatmp)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers = [36], return_contacts=True)
        token_representations = results["representations"][36]
        sequence_representations = []
        for i, (_, seq) in enumerate(datatmp):
            sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0))
        final1 = {}
        for pdb_chain,representations in zip(datatmp,sequence_representations):
            final1[pdb_chain[0]] = representations
        return final1

    complex_data = []
    for comp in seq_dict:
        for partner in seq_dict[comp]:
            for seq,seq_id in zip(seq_dict[comp][partner],range(len(seq_dict[comp][partner]))):
                complex_data.append((comp+'_'+partner+'_'+str(seq_id),seq[0]))

    for i in complex_data:
        embedding_dict = get_embed([i])
        embedding_dict[i[0]] = embedding_dict[i[0]].numpy().tolist()
        with open(f'{embedding_file_path}{i[0]}.json','w') as f:
            json.dump(embedding_dict,f)
    all_torch = {}
    for comp in seq_dict:
        comp_tensors = []
        for partner in seq_dict[comp]:
            partner_tensors = []
            for seq_id in range(len(seq_dict[comp][partner])):
                id_use = comp+'_'+partner+'_'+str(seq_id)
                with open(f'{embedding_file_path}{id_use}.json','r') as fw:
                    partner_tensor = json.load(fw)
                    partner_tensors.append(torch.tensor(list(partner_tensor.values())).squeeze(0))
            comp_tensors.append(torch.mean(torch.stack(partner_tensors), dim=0))
        all_torch[comp] = torch.cat(comp_tensors,axis = 0) 
    test_keys = [key for key in all_torch]
    test_dataset = Dataset(test_keys,all_torch)
    test_dataloader = DataLoader(test_dataset, batch_size = 32, shuffle=True) 
    label_list = {}
    node_dims = 5120
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    loss_fn = nn.MSELoss()
    df_mean = pd.DataFrame()
    for num in range(50):
        seed = 0
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        np.random.seed(seed)
        random.seed(seed) 
        torch.backends.cudnn.deterministic = True
        model_test = torch.load(pathinput + '/mlp5120_models/mlp5120_{}.pkl'.format(num),map_location=device)
        model_test.eval()
        pred_data = []
        lable_data = []
        key_data = []
        loss_data = []
        for i_batch, sample_test in enumerate(test_dataloader):
            graph = sample_test['data'].to(device)
            graph = torch.squeeze(graph , 1)
            key = sample_test['key']
            logits = model_test.forward(graph)
            pred = torch.squeeze(logits , 1)
            pred_value = torch.squeeze(logits, 1)
            for i in list(pred_value.cpu().detach().numpy()):
                pred_data.append(i)
            for i in key:
                key_data.append(i)
        tmp_df = pd.DataFrame({'key':key_data,'pred':pred_data})
        df = tmp_df.sort_values(['key']).reset_index(drop=True)
        df_mean = pd.concat([df_mean,df],axis = 1)
    df_mean['MLP_{5120}'] = np.average(df_mean['pred'], axis=1)
    df_mean['PDBid'] = df_mean['key'].iloc[:,:1]
    df_mean = df_mean[['PDBid','MLP_{5120}']]



    def model_predict():
        kinds = ['contacts','inter_sec_stru','nis','ipot','provean','McVol_whole','aaindex1','energy_asa_len']
        pdb_fea_list = []
        for kind in kinds:
            file = pdb+'_'+kind+'.'+suffix
            pdb_fea = pd.read_csv(pathoutput+'/'+file,sep='\t')
            pdb_fea_list.append(pdb_fea)
        pdb_feas = pd.concat(pdb_fea_list,axis=1)[need_feas]
        model_fit = joblib.load(path_model)
        pdb_feas['RF_{13}'] = model_fit.predict(pdb_feas)
        pdb_infos = pd.read_csv(in_file+'.cleaned',sep='\t')
        pdb_infos_feas = pd.concat([pdb_infos,pdb_feas],axis=1)
        pdb_infos_feas = pdb_infos_feas.merge(df_mean[['PDBid','MLP_{5120}']] , on = 'PDBid',how = 'left')
        pdb_infos_feas['AvgEns'] = (pdb_infos_feas['RF_{13}'] + pdb_infos_feas['MLP_{5120}']) / 2
        pdb_infos_feas.to_csv(jobpath +'/' + jobid +'.AvgEns',sep='\t',index=0)

    model_predict()
