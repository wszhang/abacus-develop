#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)
import sys
sys.stdout = Unbuffered(sys.stdout)

import os
import re
import json
import time
import copy
import math
import numpy as np
import subprocess
import string
from io import StringIO
# from string import ljust
# from _elementtree import Element
# from distutils.cygwinccompiler import get_versions
# from scipy.io.array_import import default
# #
#path_thisfile = os.path.dirname(os.path.realpath(__file__))
#sys.path.append( path_thisfile )
#sys.path.append( os.path.realpath( path_thisfile+'/../' ) )
#
##allinone avoid print buffer, instead of use sys.stdout.flush() after each print
#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
#sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)
##not work: os.environ["PYTHONUNBUFFERED"] = "1"
##print " zero print_buffer " 
##print( " ", flush=True)


def Get_FileString(file_fullPath):
    file = open( file_fullPath )
    data = file.read()
    return data

def get_string_linehead( headString, input_string ):
    linematch = re.search(r"^[ ]*"+headString+"[^#^\n]*", input_string, flags=re.DOTALL|re.MULTILINE)
    #localarray = re.split( r"[ ]+", linematch.group(), maxsplit=1 )
    localarray = linematch.group().split(maxsplit=1)
    return localarray[1].strip()

def get_nRows_linehead( headString, input_string ):
    linematch = re.findall(r"^[ ]*"+headString+"[^#^\n]*", input_string, flags=re.DOTALL|re.MULTILINE)
    #print(len(linematch))
    return len(linematch)

def get_array_linehead( headString, input_string ):
    linematch = re.search(r"^[ ]*"+headString+"[^#^\n]*", input_string, flags=re.DOTALL|re.MULTILINE)
    #localarray = re.split( r"[ ]+", linematch.group() )
    localarray = linematch.group().split()
    return localarray[1:]

def strs_to_ints( array ):
    return [ int(ii) for ii in array ]

def strs_to_floats( array ):
    return [ float(ii) for ii in array ]

def str_to_bool(str):
    return str == "True" or str == "true"

def strs_to_bools(array):
    return [ str == "True" or str == "true" for str in array ]

def Search_Num_nearStr(dftlogfpath,theStr=r'!FINAL_ETOT_IS', NumIndex=0 ):
    if os.path.isfile(dftlogfpath) :
        dftlogfile = open( dftlogfpath )
        #print "    searching num near %s in %s :"%(theStr,dftlogfpath)
        for logline in dftlogfile :
            linematch = re.search(theStr, logline)
            if linematch :
                #print logmatch.group(0) 
                print(logline.strip())
                #num_match  = re.search("[-+]?[0-9]*\.?[0-9]+",logline)
                #num_match  = re.findall("[-+]?[0-9]*\.?[0-9]+",logline)
                num_match  = re.findall(r'[-+]?\d+\.?\d*[eE]?[-+]?\d+',logline)
                #total_E = string.atof( energy_match.group(0) ) 
                #print total_E * 2 
                if num_match:
                    print(num_match)
                    #return float( num_match.group(0) ) #* 13.60569253 #13.605698
                    return float(num_match[NumIndex])
    return 0.0
#Search_Num_nearStr("./Cu/outdata/pwscf.xml",theStr=r'<etot>')
#print "%.12f"%Search_Num_nearStr("26_Fe_100/OUT.Fe-10-1.7/running_scf.log")

def parse_arguments():
    #global args
    import argparse
    # Instantiate the parser
    parser = argparse.ArgumentParser(description='Optional Description')
    
    # Required positional argument
    parser.add_argument('InputFile', type=str, nargs='?',
                    help='Absolute or relative path including the name of the input file')

    parser.add_argument('--HostList', type=str,
                    help='Host list separated by commas ')
    args = parser.parse_args()
    if args.InputFile == None:
        args.InputFile = "ORBITAL_INPUT"
    if args.HostList == None:
        args.HostList = "localhost"
    return args

def orbConf_to_list(str):
    temp_list = list( str )
    #print("temp_list:%s"%temp_list)
    results_list = []
    #print( " %20s : %s"%( "Orbitals for %2s"%element[iElement] , orbConfList_Level[iLevelm1][iElement] ), end='\n' ) 
    L=0
    for ii in list(range(0, len(temp_list), 2)): # range(len(orbConfList)/2 ):
        iip1 = ii + 1
        while temp_list[iip1] != Llabel[L]:
            results_list.append(0)
            L+=1
            # print("   ---", results_list )
        results_list.append( int(temp_list[ii]) )
        L+=1
        # print("===   ", results_list )
    # print( " Orbitals List = %s"%( results_list ), end='\n' ) 
    return results_list

def get_input_STRU( stru_type, element, mass, pseudofile, lat0, BL ):
    dis1 = BL * 0.86603 
    dis2 = BL * 0.5     
    dis3 = BL * 0.81649 
    dis4 = BL * 0.28867 
    if (stru_type == "dimer"):
        na=2
        input_STRU='''ATOMIC_SPECIES
%s %.6f %s
LATTICE_CONSTANT
%.6f  // add lattice constant(a.u.)
LATTICE_VECTORS
1 0 0
0 1 0
0 0 1
ATOMIC_POSITIONS
Cartesian_angstrom  //Cartesian or Direct coordinate.
%s      //Element Label
0.0     //starting magnetism
2       //number of atoms
0.0      0.0      0.0        0   0   0  // crystal coor.
0.0      0.0      %.6f   0   0   0
'''%(element, mass, pseudofile, lat0, element, BL)

    elif (stru_type == "trimer" ):
        na=3
        input_STRU='''ATOMIC_SPECIES
%s %.6f %s
LATTICE_CONSTANT
%.6f  // add lattice constant(a.u.)
LATTICE_VECTORS
1 0 0
0 1 0
0 0 1
ATOMIC_POSITIONS
Cartesian_angstrom  //Cartesian or Direct coordinate.
%s      //Element Label
0.0     //starting magnetism
3       //number of atoms
0.0      0.0      0.0        0   0   0  // crystal coor.
0.0      0.0      %.6f   0   0   0
0.0      %.6f %.6f   0   0   0
'''%(element, mass, pseudofile, lat0, element, BL,dis1,dis2)

    elif (stru_type == "tetramer" ): 
        na=4
        input_STRU='''ATOMIC_SPECIES
%s %.6f %s
LATTICE_CONSTANT
%.6f  // add lattice constant(a.u.)
LATTICE_VECTORS
1 0 0
0 1 0
0 0 1
ATOMIC_POSITIONS
Cartesian_angstrom  //Cartesian or Direct coordinate.
%s      //Element Label
0.0     //starting magnetism
4       //number of atoms
0.0      0.0      0.0        0   0   0  // crystal coor.
0.0      0.0      %.6f    0   0   0 
0.0      %.6f %.6f   0   0   0 
%.6f %.6f %.6f   0   0   0 
'''%(element, mass, pseudofile, lat0, element, BL,dis1,dis2,dis3,dis4,dis2)
    #print(input_STRU)
    return input_STRU, na

def get_input_KPOINTS( ):
    input_KPOINTS='''K_POINTS
0
Gamma
1 1 1 0 0 0
'''
    return input_KPOINTS

def get_input_INPUTw(STRUname, element, rcut, BL):
    input_INPUTw='''WANNIER_PARAMETERS
rcut 10
out_spillage 2
spillage_outdir OUT.%s-%s-%s-%s
'''%(element, STRUname, rcut, BL)
    return input_INPUTw

def get_input_INPUTs(ecut, rcut):
    input_INPUTs='''INPUT_ORBITAL_INFORMATION
<SPHERICAL_BESSEL>
1           // smooth or not
0.1         // smearing_sigma
%s        // energy cutoff for spherical bessel functions(Ry)
%s        // cutoff of wavefunctions(a.u.)
1.0e-12     // tolerence
</SPHERICAL_BESSEL>
'''%(ecut, rcut) 
    return input_INPUTs

def get_input_INPUT(name, Pseudo_dir, nspin, maxL, nbands_STRU, Ecut, smearing_sigma):
    input_INPUT='''INPUT_PARAMETERS
suffix              %s
stru_file           %s.stru
pseudo_dir          %s
kpoint_file         KPOINTS
wannier_card        INPUTw
calculation         scf
ntype               1
nspin               %s
lmaxmax             %s

symmetry            0
nbands             	%s

ecutwfc             %s
scf_thr             1.0e-7  // about iteration
scf_nmax            1500

smearing_method     gauss
smearing_sigma      %s

mixing_type         pulay       // about charge mixing
mixing_beta         0.4
mixing_ndim         8
printe				1
'''%(name, name, Pseudo_dir, nspin, maxL, nbands_STRU, Ecut, smearing_sigma)
    return input_INPUT

def write_string_tofile(input, filename):
    ifile = open(filename, 'w')
    ifile.write(input)
    ifile.flush()
    ifile.close()


def pw_calculation(iElement, iEcut, iRcut, STRUList):
    #iElement = 0
    #iEcut=0
    #iRcut=0
    namePW = {}

    for STRUname in STRUList:
      print( "\n %s "%( "-"*92) )
      print( " %s Get PW WaveFunction for %s with %s bond-lengths  %s"%("-"*20, STRUname, nBL_STRU[STRUname], "-"*20) )
      print( " %s "%( "-"*92) )

      if (type_STRU[STRUname] == "customer"):
        print( " Use customized PW WaveFunction Dir: " )
        for iBL in range( nBL_STRU[STRUname] ):
            print( " %60s"%datapath_STRU[STRUname][iBL] )

      elif (type_STRU[STRUname] != "customer"):
        print( "\n Setting PW WaveFunction Dir: " )
        namePW[STRUname] = [ None for iBL in range(nBL_STRU[STRUname]) ]

        for iBL in range( nBL_STRU[STRUname] ):
            namePW[STRUname][iBL] = "%s-%s-%s-%s"%(element[iElement], STRUname, Rcut[iRcut], BL_STRU[STRUname][iBL] )
            datapath_STRU[STRUname][iBL] = "../OUT.%s"%namePW[STRUname][iBL] 
            print( " %60s"%datapath_STRU[STRUname][iBL] )
        #
        for iBL in range( nBL_STRU[STRUname] ):
            print( "\n %s Do PW Calculation with Bond Length: %s %s"%('-'*25, BL_STRU[STRUname][iBL], '-'*25 ))

            (input_STRU, nAtoms) = get_input_STRU( type_STRU[STRUname], element[iElement], mass, Pseudo_name[iElement], lat0, BL_STRU[STRUname][iBL] )
            print( " %20s = %s"%("nAtoms", nAtoms), end='\n')
            # print(input_STRU, nAtoms)
            write_string_tofile(input_STRU, "%s.stru"%namePW[STRUname][iBL] )

            input_KPOINTS = get_input_KPOINTS()
            # print(input_KPOINTS)
            write_string_tofile(input_KPOINTS, "KPOINTS")

            input_INPUTw = get_input_INPUTw(STRUname, element[iElement], Rcut[iRcut], BL_STRU[STRUname][iBL])
            # print(input_INPUTw)
            write_string_tofile(input_INPUTw, "INPUTw")

            input_INPUTs = get_input_INPUTs( Ecut[iEcut], Rcut[iRcut] )
            # print(input_INPUTs)
            write_string_tofile(input_INPUTs, "INPUTs")

            input_INPUT = get_input_INPUT( namePW[STRUname][iBL], Pseudo_dir, 
                            nspin_STRU[STRUname], maxL_STRU[STRUname], nbands_STRU[STRUname], Ecut[iEcut], sigma )
            # print(input_INPUT)
            write_string_tofile(input_INPUT, "INPUT")

            PW_WF_file = "OUT.%s/orb_matrix.1.dat"%namePW[STRUname][iBL]
            sys_run_str = '''%s;
pwd;
which mpirun mpiexec.hydra;
if [ ! -f "%s" ]; then
    echo " Not found %s, Calculating PW WF ... "
    %s %s;
else
    echo " Has found %s, Skip Calculation "
fi
#mpiexec.hydra -np 4 -host localhost /gpfs/home/nic/wszhang/abacus222_intel-2018u4/ABACUS.mpi
'''%(EXE_bash_env,  PW_WF_file, PW_WF_file, EXE_mpi, EXE_pw, PW_WF_file)
            print(" runcmd: \n%s \n"%sys_run_str )
            sys.stdout.flush() 

            # os.environ["PYTHONUNBUFFERED"] = "1"
            #process = subprocess.Popen(["your_cmd"]...)
            #process.wait()
            #with open('env_ff.txt', 'w') as ff:
            #    out = subprocess.run('env', shell=True, stdout=ff, text=True)
            #subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL, timeout=7200) 

            subprocess.run( [sys_run_str, "--login"], shell=True, text=True, timeout=72000) 
            sys.stdout.flush() 

            #osvalue = os.system(sys_run_str) 
            #print(" osvalue = %s"%osvalue )
            #if osvalue != 0 : print(" fail to call "+sys_run_str) #continue #quit() 
            #sys.stdout.flush() 
            #break 


##################################  Prepare SIAB INPUT ##################################
def Prepare_SIAB_INPUT(iEcut, iRcut, iLevel):
    iLevelm1 = iLevel-1
 
    STRUname = refSTRU_Level[iLevelm1]
    print( "\n %s "%( "-"*92) )
    print( " %s Do SIAB Calculation for Level%s with Ref %-10s %s"%("-"*20, iLevel, STRUname, "-"*20) )
    print( " %s "%( "-"*92) )

    INPUT_json = {"file_list":{}, "info":{}, "weight":{}, "C_init_info":{}, "V_info": {} }

    INPUT_json["file_list"] = {"origin":[], "linear":[] }
    INPUT_json["file_list"]["origin"] = [ datapath_STRU[STRUname][iBL]+"/orb_matrix.0.dat" for iBL in range(nBL_STRU[STRUname]) ]
    INPUT_json["file_list"]["linear"] = [ [ datapath_STRU[STRUname][iBL]+"/orb_matrix.1.dat" for iBL in range(nBL_STRU[STRUname]) ] ]

    INPUT_json["info"] = {"Nt_all": element, 
			"Nu":   { element[iElement]:orbConf_to_list(orbConf_Level[iLevelm1][iElement]) for iElement in range(len(element) )  },
			"Rcut": { element[iElement]:Rcut[iRcut] for iElement in range(len(element)) },
			"dr":   { element[iElement]:0.01 for iElement in range(len(element)) },
			"Ecut": { element[iElement]:int(Ecut[iEcut]) for iElement in range(len(element)) }, 
            "lr": 0.01, 
            "cal_T": True,  "cal_smooth": True, "max_steps": max_steps } 
    
    INPUT_json["weight"] = { "stru": [1] * nBL_STRU[STRUname], 
                             "bands_range": refBandsRange_Level[iLevelm1] }
    
    INPUT_json["C_init_info"] = {}
    if ( fixPre_Level[iLevelm1] == "None" or fixPre_Level[iLevelm1] == "none" ):
        INPUT_json["C_init_info"]["init_from_file"] = False
    else:
        INPUT_json["C_init_info"]["init_from_file"] = True

        if  (iLevelm1  > 0):
            INPUT_json["C_init_info"]["C_init_file"] = "Level%s.ORBITAL_RESULTS.txt"%iLevelm1
        elif(iLevelm1 == 0):
            INPUT_json["C_init_info"]["C_init_file"] = "ORBITAL_RESULTS.txt"

        if ( fixPre_Level[iLevelm1] == "fix" or fixPre_Level[iLevelm1] == "Fix"  ):
            INPUT_json["C_init_info"]["opt_C_read"] = False
        else:
            INPUT_json["C_init_info"]["opt_C_read"] = True

    INPUT_json["V_info"] = {
        "same_band":True,
        "init_from_file":True }

    if ( 'opt_C_read' in INPUT_json["C_init_info"] ):
        if ( INPUT_json["C_init_info"]["opt_C_read"] == True  and  INPUT_json["C_init_info"]["init_from_file"] == True ):
            INPUT_json["info"]["lr"] = 0.001

    INPUT_json_str = json.dumps(INPUT_json, indent=2)
    #print(INPUT_json_str)

    ifile_input = open("INPUT", 'w')
    ifile_input.write(INPUT_json_str)
    ifile_input.flush()
    ifile_input.close()


###################################  Setting Constants ###################################
Hartree_to_eV=27.21138505
periodtable = {   'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
                  'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
               'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
               'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
               'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
               'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
               'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
               'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
               'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
               'Ba': 56, #'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
                    ## 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
                    ## 'Er': 68, 'Tm': 69, 'Yb': 70, 
                    ## 'Lu': 71, 
               'Hf': 72, 'Ta': 73,
               'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
               'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 
                    ## 'Po': 84, #'At': 85,
                    ## 'Rn': 86, #'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
                    ## 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
                    ## 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
                    ## 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
                    ## 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Uut': 113,
                    ## 'Fl': 114, 'Uup': 115, 'Lv': 116, 'Uus': 117, 'Uuo': 118
               } 
periodtable_r = {v: k for k, v in periodtable.items()}
#print len(periodtable_r.keys() )
#print periodtable_r[1]
#print periodtable['Hg']
Llabel=['s','p','d','f','g']
lat0=20.0
mass=1.0


###################################  Parse Arguments ###################################
args=parse_arguments()
print(" Starting SIAB\n")
print( " %s "%('*'*92) )
print( " *%s*"%(" "*90 ) )
print( " * %75s %s *"%("Systematically Improvable Atomic-orbital Basis (SIAB) generator", " "*12 )  )
print( " * %65s %s *"%("for Linear Combination of Atomic Orbitals (LCAO)", " "*22 )  )
print( " *%s*"%(" "*90) )
print( " %s "%('*'*92) )
print("\n %20s = %s "%("InputFile",args.InputFile))
#print(" HostList: %s "%args.HostList)


###################################  Parse InputFile ###################################
input_string=Get_FileString(args.InputFile)
#print(input_string)
EXE_bash_env= get_string_linehead( "EXE_bash_env", input_string )
EXE_mpi     = get_string_linehead( "EXE_mpi", input_string )
EXE_pw      = get_string_linehead( "EXE_pw", input_string )
EXE_orbital = get_string_linehead( "EXE_orbital", input_string )
element     = get_array_linehead( "element", input_string )
Ecut        = strs_to_ints(get_array_linehead( "Ecut", input_string ) )
Rcut        = strs_to_ints(get_array_linehead( "Rcut", input_string ) )
sigma       = float(get_string_linehead( "sigma", input_string ) )
Pseudo_dir  = get_string_linehead( "Pseudo_dir", input_string )
Pseudo_name = get_array_linehead( "Pseudo_name", input_string )
max_steps   = int(get_string_linehead( "max_steps", input_string ) )

element_num = [ periodtable[ii] for ii in element ]

print(" %20s = %s "%("EXE_bash_env", EXE_bash_env) )
print(" %20s = %s "%("EXE_mpi", EXE_mpi) )
print(" %20s = %s "%("EXE_pw", EXE_pw) )
print(" %20s = %s "%("EXE_orbital", EXE_orbital) )
print(" %20s = %s "%("element", element) )
print(" %20s = %s "%("element_num", element_num) )
print(" %20s = %s "%("Ecut", Ecut) )
print(" %20s = %s "%("Rcut", Rcut) )
print(" %20s = %s "%("sigma", sigma) )
print(" %20s = %s "%("Pseudo_dir", Pseudo_dir) )
print(" %20s = %s "%("Pseudo_name", Pseudo_name) )
print(" %20s = %s "%("max_steps", max_steps) )

input={}

nLevel = get_nRows_linehead( "Level", input_string )
print("\n %20s = %s "%("nLevel", nLevel) )
for iLevel in range(1,nLevel+1):
    input["Level%s"%iLevel] = get_array_linehead( "Level%s"%iLevel, input_string )
    #print(" %20s : %s"%( "Level%s"%iLevel, input["Level%s"%iLevel] ) )

refSTRU_Level=[]
fixPre_Level=[]
maxL_Level=[]
orbConf_Level=[]
restartLevel=[]
refBands_Level=[]

for iLevel in range(1,nLevel+1):
    iLevelm1 = iLevel - 1

    #print( "\n Level%s: "%iLevel, end='\n ')
    print( " %-6s:"%("Level%s"%iLevel), end='\n')

    refSTRU_Level.append( input["Level%s"%iLevel][0] )
    print( " %20s = %s"%("Reference Struture", refSTRU_Level[iLevelm1] ), end='\n')

    refBands = eval(input["Level%s"%iLevel][1])
    refBands_Level.append( refBands )
    print( " %20s = %s"%("Reference Bands", refBands_Level[iLevelm1]), end='\n')

    fixPre_Level.append( input["Level%s"%iLevel][2] )
    print( " %20s : %s"%("Fix input orbitals?", fixPre_Level[iLevelm1]), end='\n')

    if iLevel == 1 :
        restartLevel.append(False)
    else:
        restartLevel.append(True)
    print( " %20s : %s"%("Restart Level?", restartLevel[iLevelm1]), end='\n') 

    orbConf_Level.append( input["Level%s"%iLevel][3:] )
    # print( " %20s = %s"%("Orbital Conf", orbConf_Level[iLevelm1]), end='\n') 
    for iElement in range(len(element) ):
        print( " %20s : %s / %s"%( "Orbitals for %2s"%element[iElement] , \
                    orbConf_Level[iLevelm1][iElement], \
                    orbConf_to_list(orbConf_Level[iLevelm1][iElement]) ), end='\n' ) 

    print( " ", end='\n' ) 
# orbConfList_Level = copy.deepcopy(orbConf_Level)
# for iLevel in range(1,nLevel+1):
#     iLevelm1 = iLevel - 1 
#     for iElement in range(len(element) ):
#         orbConfList_Level[iLevelm1][iElement] = orbConf_to_list(orbConf_Level[iLevelm1][iElement])
#         print( " %20s : %s"%( "Orbitals for %2s"%element[iElement] , orbConfList_Level[iLevelm1][iElement] ), end='\n' ) 


STRUList = list(dict.fromkeys(refSTRU_Level))
nSTRU = len(STRUList)
#nSTRU = get_nRows_linehead( "STRU", input_string )
print("\n Parse %s types of structures: %s"%(nSTRU, STRUList) )

type_STRU={}
maxL_STRU={}
nspin_STRU={}
BL_STRU={}
nBL_STRU={}
datapath_STRU={}
nbands_STRU={}

for STRUname in STRUList:
    input[STRUname] = get_array_linehead( STRUname, input_string )
    #print(" %20s : %s"%( STRUname, input[STRUname] ) )
    print( " %s:"%(STRUname), end='\n')

    type_STRU[STRUname] = input[STRUname][0]
    print( " %20s = %s"%("STRU Type", type_STRU[STRUname]), end='\n')

    if (type_STRU[STRUname] != "customer"):
        nbands_STRU[STRUname] = int(input[STRUname][1])
        print( " %20s = %s"%("nbands", nbands_STRU[STRUname]), end='\n')

        maxL_STRU[STRUname] = int(input[STRUname][2])
        print( " %20s = %s"%("maxL", maxL_STRU[STRUname] ), end='\n')

        nspin_STRU[STRUname] = int(input[STRUname][3])
        print( " %20s = %s"%("nspin", nspin_STRU[STRUname]), end='\n')

        BL_STRU[STRUname] = strs_to_floats(input[STRUname][4:]) 
        print( " %20s = %s"%("Bond Length List", BL_STRU[STRUname]), end='\n')

        nBL_STRU[STRUname] = len(BL_STRU[STRUname])
        #print(  " %20s = %s"%("BL List Size", nBL_STRU[STRUname]), end='\n')

        datapath_STRU[STRUname] = ["None" for ii in range(nBL_STRU[STRUname]) ]
    else:
        datapath_STRU[STRUname] = input[STRUname][1:]
        print(  " %20s = %s"%("WF Data Path", datapath_STRU[STRUname]), end='\n')

        nBL_STRU[STRUname] = len(datapath_STRU[STRUname])
    
    print(  " %20s = %s"%("BL List Size", nBL_STRU[STRUname]), end='\n')

nSave = get_nRows_linehead( "Save", input_string )
if (nSave > 0):
    print("\n %20s : %s "%("Save Orbital?", True if nSave>0 else False) )
    print(  " %20s = %s "%("nSave", nSave) )
    for ii in range(1,nSave+1):
        input["Save%s"%ii] = get_array_linehead( "Save%s"%ii, input_string )
        print(" %s %s as %s"%( "Save", input["Save%s"%ii][0] , input["Save%s"%ii][1]) )
print( " ", end='\n' ) 


##################################  Derived parameter  ##################################
nRcut = len(Rcut)

refBandsRange_Level=[]
for iLevel in range(1,nLevel+1):
    iLevelm1 = iLevel - 1
    STRUname = refSTRU_Level[iLevelm1]
    #print( type(refBands_Level[iLevelm1] ) )

    if ( type(refBands_Level[iLevelm1]) == list ):
        refBandsRange_Level.append( refBands_Level[iLevelm1] )
    else:
        refBandsRange_Level.append( [ refBands_Level[iLevelm1] ] * nBL_STRU[STRUname] )
    #print( iLevelm1, refBandsRange_Level[iLevelm1])


##################################  Do    Calculation ##################################
iElement=0
iEcut=0

for iRcut in range(nRcut):
    ################################  Do PW Calculation ################################
    ElementDir = os.getcwd()
    print(" Current working directory %s "%ElementDir )
    # Now change the directory
    SIAB_wdir = "%s_%s_%sRy"%(element_num[iElement], element[iElement], Ecut[iEcut])
    try:
        os.mkdir(ElementDir+"/"+SIAB_wdir)
    except OSError as error:
        print(" Already has directory: %s"%(ElementDir+"/"+SIAB_wdir) )    
    SIAB_fullpath= ElementDir+"/"+SIAB_wdir
    os.chdir(SIAB_fullpath)
    # Check current working directory.
    print(" Directory changed to %s "%os.getcwd() )
    pw_calculation(iElement, iEcut, iRcut, STRUList)


    ################################  Do SIAB Calculation ###############################
    for iLevel in range(1,nLevel+1):
        iLevelm1 = iLevel - 1 

        print(" Current working directory %s "%os.getcwd() )
        SIAB_rcutdir = SIAB_fullpath+"/%s"%(Rcut[iRcut])
        try:
            os.mkdir(SIAB_rcutdir)
        except OSError as error:
            print(" Already has directory: %s"%SIAB_rcutdir )    
        os.chdir(SIAB_rcutdir)
        # Check current working directory.
        print(" Directory changed to %s "%os.getcwd() )
        Prepare_SIAB_INPUT(iEcut, iRcut, iLevel)

        sys_run_str = '''
#conda init bash

%s;
nproc_pw=`%s hostname| wc -l`
export OMP_NUM_THREADS=$nproc_pw
echo " OMP_NUM_THREADS:" $OMP_NUM_THREADS

#module purge
#module load anaconda3_nompi
#module load gcc/9.2.0
#module list 2>&1;
pwd;
#source activate pytorch110 
#conda activate pytorch110 
#conda info --envs
echo python: `which python3`

python3 %s

#conda deactivate
'''%(EXE_bash_env, EXE_mpi, EXE_orbital)

        subprocess.run( [ sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL, timeout=18000) 
        sys.stdout.flush() 

        Leveln = "Level"+str(iLevel)
        sys_run_str = '''
mv INPUT                %s.INPUT
mv ORBITAL_%sU.dat      %s.ORBITAL_%sU.dat
mv ORBITAL_PLOTU.dat    %s.ORBITAL_PLOTU.dat
mv ORBITAL_RESULTS.txt  %s.ORBITAL_RESULTS.txt
mv Spillage.dat         %s.Spillage.dat
'''%( Leveln, element_num[iElement],  Leveln,element_num[iElement],  Leveln,   Leveln,  Leveln ) 
        subprocess.run( [ sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL, timeout=60) 


    ##################################  Save Orbitals  #################################
    print( "\n %s "%( "-"*92) )
    print( " %s           %s           %s"%("-"*26, " Save Orbitals ", "-"*26) )
    print( " %s "%( "-"*92) )

    print("\n Current working directory %s "%os.getcwd() )
    # Now change the directory
    os.chdir(ElementDir)
    # Check current working directory.
    print(" Directory changed to %s "%os.getcwd() )

    if (nSave > 0):
        for ii in range(1,nSave+1):
            Leveln = input["Save%s"%ii][0]
            Leveln = int(Leveln[5:])
            resultPath = SIAB_wdir+"/"+str(Rcut[iRcut])+"/Level"+str(Leveln) 
            orbName = element[iElement]+"_gga_"+str(Rcut[iRcut])+"au_"+str(Ecut[iEcut])+"Ry_"+orbConf_Level[Leveln-1][iElement]

    # need to check the orbConf_Level in files
    # Number of Sorbital-->       3
    # Number of Porbital-->       3
    # Number of Dorbital-->       2

            orbType = input["Save%s"%ii][1]
            try:
                os.mkdir(orbType)
            except OSError as error:
                print(" Already has directory: %s"%( orbType ) ) 
            savePath   = orbType+"/"+str(Rcut[iRcut])
            try:
                os.mkdir(savePath)
            except OSError as error:
                print(" Already has directory: %s"%( savePath ) ) 
            print("\n Save Level%s results to dir: %s"%(str(Leveln), savePath) )

            sys_run_str = '''
cp -avp %s.INPUT                %s/INPUT
cp -avp %s.ORBITAL_%sU.dat      %s/%s.orb
cp -avp %s.ORBITAL_PLOTU.dat    %s/ORBITAL_PLOTU.dat
cp -avp %s.ORBITAL_RESULTS.txt  %s/ORBITAL_RESULTS.txt
cp -avp %s.Spillage.dat         %s/Spillage.dat
'''%( resultPath, savePath, resultPath,element_num[iElement],orbType,orbName, resultPath,savePath, resultPath,savePath, resultPath,savePath ) 

            subprocess.run( [ sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL, timeout=60) 

        print("\n Saved %s Orbitals"%(nSave) )
    print( " ", end='\n' ) 
