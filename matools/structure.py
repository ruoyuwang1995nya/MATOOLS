# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:00:41 2022

@author: 王若宇
"""

import numpy as np
#import sys
import yaml
from matools import vector,utility, atom



# ======================================================        
        
class ionsClassIter():
    '''
    Class iterator for "ions" cls
    '''
    def __init__(self,ions):
        self._ions=ions # take ions instance as input
        self._n=0     # Current index for iteration
        
        
    def __iter__(self):
        '''Return iterator object'''
        return self
    
    def __next__(self):
        if self._n < len(self._ions._ions):
            member=self._ions._ions[self._n]
            self._n += 1
            return member
        else:
            raise StopIteration
        
class ions_cell():
    '''
    Class for managing the group of ions within unit cell
    '''
    def __init__(self):
        self._ions=[]  # include atom and position within cell
        
    def __getitem__(self,index:int): # make ions_cell object subscriptable
        return self._ions[index]
        
    def members(self):
        return [(index,i[0].element,i[1].vec) for index,i in enumerate(self._ions)]
        
    def compact(self):
        '''
        return a compact list of elements and atom number
        '''
        if len(self._ions)!=0:
            last_ele=self._ions[0][0].element
            last_num=1
            compact_ls=[]
            for i in range(1,len(self._ions)):
                if self._ions[i][0].element==last_ele:
                    last_num+=1
                else:
                    compact_ls.append((last_ele,last_num))
                    last_ele=self._ions[i][0].element
                    last_num=1
            compact_ls.append((last_ele,last_num))
            return compact_ls
        else:
            print("No atoms")
          
            
        
        
    def add_ion(self,ion:atom.atom,position:vector.lat_vec):
        self._ions.append((ion,position))
        
    def remove_ion(self,index:int):
        '''
        Remove designated atom from ions group
        '''
        if index < len(self._ions):
            ion=self._ions[index]
            #print()
            self._ions.remove(ion)
            print('{} at {} Removed'.format(ion[0].element,ion[1].vec))
        else:
            raise IndexError
        
    def __iter__(self):
        '''Return iterator object'''
        return ionsClassIter(self)
        
    def __len__(self):
        '''return the length of object'''
        return len(self._ions)
        
    def generator(self):
        '''generate ions on the flow'''
        i=0
        while i < len(self):
            yield self[i]
            i+=1
    
## ====================== structure ============================
class structure():
    '''
    Class for crystal structure
    
    Unit length: Angstrom    
    '''
    def __init__(self, name:str, a:vector.lat_vec,b:vector.lat_vec,c:vector.lat_vec, ions, coord_sys,scaler):
            #print("empty structure")
        self.name=name
        self.scaler=scaler

        # set the lattice vectors in Cartesian coordinate
        self.a=a      
        self.b=b
        self.c=c
                
        self._coord_sys=coord_sys
        self._ions=ions

            
    @property
    def coord_sys(self):
        return self._coord_sys
    @coord_sys.setter
    def coord_sys(self,sys:str):
        if sys[0] not in 'CDcd':
            raise TypeError("Invalid coordinate system, must be cartesian or direct system")
        else:
            self._coord_sys=sys
            
    @property
    def ions(self):
        return self._ions
    @ions.setter
    def ions(self,value:ions_cell):
        self._ions=value

# ======================== I/O interface =========================
    # ==================== Input file =========================
    @classmethod
    def from_poscar(cls,filename):
        '''read cell structure from POSCAR'''
        try:
            with open(filename,'r',encoding='utf-8') as f:
                lines=f.readlines()

            name=lines[0].strip('\n')
            scaler=float(lines[1].strip('\n'))
            a=vector.lat_vec(np.array(lines[2].strip('\n').split()).astype('float64'))
            b=vector.lat_vec(np.array(lines[3].strip('\n').split()).astype('float64'))
            c=vector.lat_vec(np.array(lines[4].strip('\n').split()).astype('float64'))
            ele_keys=np.array(lines[5].strip('\n').split())
            ele_num=np.array(lines[6].strip('\n').split()).astype('int')
            
            ions=[]
            for atom_ele in zip(ele_keys,ele_num):
                for i in range(atom_ele[1]):
                    ions.append(atom.atom(atom_ele[0]))
            ions=np.array(ions)
            
            coordination_sys=lines[7].strip('\n')
            
            # ion positions
            ions_pos=[]
            ion_position=lines[8:8+ele_num.sum()]
            #ions_pos=np.zeros((len(ion_position),3))
            for i,n in enumerate(ion_position):
                line=n.split()
                ions_pos.append(vector.lat_vec([float(line[j]) for j in range(3)]))

            ions_cell_c=ions_cell()
            for ion in zip(ions,ions_pos):
                ions_cell_c.add_ion(ion[0],ion[1])
            
            return cls(name, a, b, c, ions_cell_c,coordination_sys,scaler=scaler)
        except OSError:
            print('cannot open', filename)

    @classmethod     
    def from_qe(cls,filename):
        pass

    @classmethod
    def from_phonopy(cls,filename):
        try:
            with open(filename,'r',encoding='utf-8') as f:
                output=yaml.full_load(f)
            name='default'
            scaler=1.0
            a=vector.lat_vec(np.array(output['lattice'][0]).astype('float64'))
            b=vector.lat_vec(np.array(output['lattice'][1]).astype('float64'))
            c=vector.lat_vec(np.array(output['lattice'][2]).astype('float64'))
            points=output['points']
            #num_pts=len(points)
            ions=[]
            ions_pos=[]
            for ion in points:
                ions.append(atom.atom(ion['symbol']))
                ions_pos.append(vector.lat_vec(ion['coordinates']))
            
            ions_cell_c=ions_cell()
            for ion in zip(ions,ions_pos):
                ions_cell_c.add_ion(ion[0],ion[1])
            
            coordination_sys='Direct'
            print(scaler)
            return cls(name, a, b, c, ions_cell_c,coordination_sys,scaler=scaler)
            
        except OSError:
            print('cannot open', filename)
            
    # ========================== Output file ======================  
    def write_to_poscar(self, output_f, name='default'):
        '''write_to_POSCAR'''
        with open(output_f,'w') as output_f:
            output_f.write(name+'\n'+str(self.scaler)+'\n')
            output_f.write(utility.iterative_print(self.a.vec)+'\n')
            output_f.write(utility.iterative_print(self.b.vec)+'\n')
            output_f.write(utility.iterative_print(self.c.vec)+'\n')

            ions=self.ions.members()
            
            ele_key=[ions[0][1]]
            ele_num=[]
            c_ele_num=1 # number of current element
            for i in range(1,len(ions)):
                if ions[i][1]==ele_key[-1]:
                    c_ele_num+=1
                else:
                    ele_key.append(ions[i][1])
                    ele_num.append(c_ele_num)
                    c_ele_num=1
            ele_num.append(c_ele_num)

            output_f.write(utility.iterative_print(ele_key)+'\n')
            output_f.write(utility.iterative_print(ele_num)+'\n')
            output_f.write(self.coord_sys+'\n')
            
            # adding element counter
            for i in ions:                               # atom position
                output_f.write(utility.iterative_print(i[2])+'\n')
        return output_f      

# ========================== Functionalities =======================
    def volume(self):
        '''
        return the volume of cell in Angstrom**3
        
        a X b . c
        '''
        a=self.a.vec
        b=self.b.vec
        c=self.c.vec
        return (a[1]*b[2]-a[2]*b[1])*c[0]+(a[2]*b[0]-a[0]*b[2])*c[1]+(a[0]*b[1]-a[1]*b[0])*c[2]
    
    def coord_transform(self):
        '''
        Transform the coordinates from Direct to Cartesian or vice versa
        '''
        if self.coord_sys[0] in 'Cc':
            for ion in self.ions:
                self.cartesian_to_direct(ion[1])
            self.coord_sys='Direct'
            print("Changed to Direct coordinates")
            
        else:
            for ion in self.ions:
                self.direct_to_cartesian(ion[1])
            self.coord_sys='Cartesian'
            print("Changed to Cartesian coordinates")
        
    def cartesian_to_direct(self,ion:vector.lat_vec):
        a=vector.multiply(ion,self.a)/self.a.modulus()**2
        b=vector.multiply(ion,self.b)/self.b.modulus()**2
        c=vector.multiply(ion,self.c)/self.c.modulus()**2
        ion.vec[0]=a
        ion.vec[1]=b
        ion.vec[2]=c
        #print("Changed to direct coordinates")
    
    def direct_to_cartesian(self,ion:vector.lat_vec):
        '''
        transform vector in direct sys to cartesian 
        '''
        #sum([])
        x=ion.vec[0]*self.a.vec[0]+ion.vec[1]*self.b.vec[0]+ion.vec[2]*self.c.vec[0]
        y=ion.vec[0]*self.a.vec[1]+ion.vec[1]*self.b.vec[1]+ion.vec[2]*self.c.vec[1]
        z=ion.vec[0]*self.a.vec[2]+ion.vec[1]*self.b.vec[2]+ion.vec[2]*self.c.vec[2]
        ion.vec[0]=x
        ion.vec[1]=y
        ion.vec[2]=z
        #print("Changed to cartesian coordinates")
            
    
    def bond_length(self,index1:int,index2:int):
        '''
        return the bond length between two atoms
        param: index
        '''
        vec=self.ions[index1][1]-self.ions[index2][1]
        if self.coord_sys[0] in 'Cc':
            return vec.modulus()
        else:
            self.direct_to_cartesian(vec)
            return vec.modulus()
            
    def nearest_neighor(self,atom:int):
        '''
        return the nearest neighor of selected atom
        '''
        try:
            pos=0
            if atom==0:
                distance=self.bond_length(0,1)
            else:
                distance=self.bond_length(0,atom)
            for i in range(len(self.ions)):
                if i == atom:
                    pass
                else:
                    distance_tmp=self.bond_length(atom,i)
                    if distance > distance_tmp:
                        distance=distance_tmp
                        pos=i
            return self.ions.members()[pos],distance
        except IndexError:
            print("atom out of index")
            
    def neighors(self,atom:int,tolerance=0.1):
        '''
        search for the nearest neigbors of selected atom within given tolerance
        '''
        try:
            #neighor_atom=[]
            distance=self.nearest_neighor(atom)[1]
            return [self.ions.members()[i] for i in range(len(self.ions)) if atom != i and self.bond_length(atom,i)< (distance+tolerance)]
        
        except IndexError:
            print("atom out of index")
            
    
    def bond_search (self,atom:int,max_length=1):
        '''
        Search atoms within max_length (in angstrom) of the central atom
        param:
        atom: index of the central atom
        max_length: the cut-off distance, in angstrom
        '''
        try:
            return [self.ions.members()[i] for i in range(len(self.ions)) if atom != i and self.bond_length(atom,i)< max_length
                    
        except IndexError:
            print("atom out of index")
    

    def bond_angle(self,atom_center:int,atom1:int,atom2:int):
        '''
        Return bond angle between two atoms adjacent to a central atom
        wait to be finished  
        '''
        pass
        
