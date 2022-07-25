import numpy as np
import sys
import yaml
import matools.vector as vector
        
## ====================== structure ============================
            
class structure():
    '''
    Class for cell structure
    Unit length: Angstrom    
    '''
    def __init__(self, name:str, a:vector.lat_vec,b:vector.lat_vec,c:vector.lat_vec, ions, ions_pos, coord_sys,scaler):
            #print("empty structure")
        self.name=name
        self.scaler=scaler

        # set the lattice vectors in Cartesian coordinate
        self.a=a      
        self.b=b
        self.c=c
                
        self._coord_sys=coord_sys
        self.ions=ions
        self.ions_pos=ions_pos

            
    @property
    def coord_sys(self):
        return self._coord_sys
    @coord_sys.setter
    def coord_sys(self,sys:str):
        if sys[0] not in 'CDcd':
            raise TypeError("Invalid coordinate system, must be cartesian or direct system")
        else:
            self._coord_sys=sys
            
        
# ======================== POSCAR generation =========================
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
            for atom in zip(ele_keys,ele_num):
                for i in range(atom[1]):
                    ions.append(atom[0])
            ions=np.array(ions)
            
            coordination_sys=lines[7].strip('\n')
            
            # ion positions
            ion_position=lines[8:8+ele_num.sum()]
            ions_pos=np.zeros((len(ion_position),3))
            for i,n in enumerate(ion_position):
                line=n.split()
                for j in range(3):
                    ions_pos[i,j]=line[j]
            ions_pos=ions_pos.astype('float64')
            return cls(name, a, b, c, ions,ions_pos,coordination_sys,scaler=scaler)
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
                ions.append(ion['symbol'])
                ions_pos.append(ion['coordinates'])
            ions=np.array(ions)
            ions_pos=np.array(ions_pos)
            coordination_sys='Direct'
            print(scaler)
            return cls(name, a, b, c, ions,ions_pos,coordination_sys,scaler=scaler)
            
        except OSError:
            print('cannot open', filename)

    
# ========================== POSCAR processing =======================
    def atom(self,n):
        '''
        return the element and position of the nth atom
        '''
        try:
            return (self.ions[n-1],self.ions_pos[n-1])
        except IndexError:
            print('atom index out of range')
    
    def bond_length(self,atom1,atom2):
        '''
        return the bond length between two atoms
        '''
        try:
            a=(self.ion[atom1-1],self.ion_position[atom1-1])
            b=(self.ion[atom2-1],self.ion_position[atom2-1])
            diff=a[1]-b[1]
            return sum(diff**2)**0.5
        except IndexError:
            print('atom index out of range')

    
    def bond_angle(atom1,atom2,atom3):
        pass
    
    def bond_search(atom,max_length):
        pass
    
    
    def iterative_print (self,ite_obj):
        '''print iteratives as string'''
        output=''
        for i in ite_obj:
            if type(i)==np.float64:
                output += (r' {:.16f}'.format(i))
            else:
                output += (r' {}'.format(i))
        return output

# ==========================   output ======================================    
    def write_to_poscar(self, output_f, name='default'):
        '''write_to_POSCAR'''
        with open(output_f,'w') as output_f:
            output_f.write(name+'\n'+str(self.scaler)+'\n')
            output_f.write(self.iterative_print(self.a)+'\n')
            output_f.write(self.iterative_print(self.b)+'\n')
            output_f.write(self.iterative_print(self.c)+'\n')
            
            num_ions=len(self.ions)
            ele_key=[self.ions[0]]
            ele_num=[]
            c_ele_num=1
            for i in range(1,len(self.ions)):
                if self.ions[i]==ele_key[-1]:
                    c_ele_num+=1
                else:
                    ele_key.append(self.ions[i])
                    ele_num.append(c_ele_num)
                    c_ele_num=1
            ele_num.append(c_ele_num)

            output_f.write(self.iterative_print(ele_key)+'\n')
            output_f.write(self.iterative_print(ele_num)+'\n')
            output_f.write(self.coordination_sys+'\n')
            for i in self.ions_pos:
                output_f.write(self.iterative_print(i)+'\n')
        return output_f
        


