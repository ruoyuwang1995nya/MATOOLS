#!/usr/bin/env python
import matools as mat
import matools.structure as structure
import numpy as np
import sys
import argparse

parser=argparse.ArgumentParser(description='A simple script generating CONTROL from POSCAR')
parser.add_argument("-ng","--ngrid",nargs='+',type=int,default=[10,10,10],help="The size of k-grid mesh")
parser.add_argument("-s","--scell",nargs='+',type=int,default=[1,1,1],help="The size of supercell used for 2FC calculation")
parser.add_argument("-t","--temperature",nargs='+',type=int,default=[300,100,400],help="Tmin, Tstep and Tmax")
parser.add_argument('--isotope',dest='isotope',action='store_true')
parser.add_argument('--autoisotope',dest='autoisotpe',action='store_true')
parser.set_defaults(isotope=False)
parser.set_defaults(autoisotope=False)
args=parser.parse_args()

if len(args.scell)==3:
    pass
else:
    print("invalid supercell dimension")
    exit(1)


input_poscar=r'./POSCAR'
output=r'./CONTROL'
#poscar=p.poscar()
try:
    poscar=structure.structure.from_poscar(input_poscar)
except:
    print("The POSCAR is invalid!")
    exit(1)

with open(output,'w') as output_f:
    ## allocation part
    output_f.write(r"&allocations"+'\n')
    output_f.write("nelements={0}\n".format(len(poscar.ions.compact())))
    output_f.write("natoms={0}\n".format(len(poscar.ions)))
    output_f.write("ngrid(:)={} {} {}\n".format(args.ngrid[0],args.ngrid[1],args.ngrid[2]))
    output_f.write(r"&end"+'\n')
    output_f.write('\n')
    ## crystal part
    output_f.write(r"&crystal"+'\n')
    # unit: angstrom
    output_f.write("lfactor = 0.1\n")
    # write the lattice vector
    output_f.write("lattvec(:,1)="+mat.utility.iterative_print(poscar.a.vec)+',\n')
    output_f.write("lattvec(:,2)="+mat.utility.iterative_print(poscar.b.vec)+',\n')
    output_f.write("lattvec(:,3)="+mat.utility.iterative_print(poscar.c.vec)+',\n')
    
    # write the elements and numbers
    ions=poscar.ions
    ele_number=ions.compact()
    elements=[]
    types=[]
    for i in ele_number:
        elements.append(i[0])
        types.append(i[1])
        
    # writing
    output_f.write("elements="+mat.utility.iterative_print(elements)+',\n')
    output_f.write("types="+mat.utility.iterative_print(types)+',\n')
    
    # write the atom positions
    for num, i in enumerate(ions):
        #position_str=poscar.iterative_print(position)
        output_f.write("positions(:,{0}) = {1}".format(num+1,mat.utility.iterative_print(i[1].vec))+'\n')
        
    # Write the supercell dimension
    output_f.write("scell(:)= ")
    for i in args.scell:
        output_f.write("{} ".format(i))
    output_f.write('\n')
    output_f.write(r"&end"+'\n')
    output_f.write('\n')
    
    ## parameters part
    output_f.write(r"&parameters"+'\n')
    output_f.write("T_min={} \n".format(args.temperature[0]))
    output_f.write("T_max={} \n".format(args.temperature[2]))
    output_f.write("T_step={} \n".format(args.temperature[1]))
    output_f.write(r"&end"+'\n')
    output_f.write('\n')
    ## flags part
    
    output_f.write(r"&flags"+'\n')
    if args.isotope:
        output_f.write("isotopes=.true.\n")
        if not args.autoisotope:
            output_f.write("autoisotopes=.false.\n")
    else:
        output_f.write("isotopes= .false.\n")
    output_f.write("convergence= .true.\n")
    output_f.write(r"&end"+'\n')

print('done')
