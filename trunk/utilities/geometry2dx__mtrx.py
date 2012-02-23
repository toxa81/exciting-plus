#!/usr/bin/python

#
# GEOMTERY.OUT to OpenDX 
#  
# Created: April 2009 (AVK)
# Modified: February 2012 (AVK) 
#

import math
import sys

def r3minv(a, b):
    t1 = a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]) + \
         a[0][1] * (a[1][2] * a[2][0] - a[1][0] * a[2][2]) + \
         a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
    if math.fabs(t1) < 1e-40:
        print "r3mv: singular matrix"
        sys.exit(0)
    
    t1 = 1.0/t1
    
    b[0][0] = t1 * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
    b[0][1] = t1 * (a[0][2] * a[2][1] - a[0][1] * a[2][2])
    b[0][2] = t1 * (a[0][1] * a[1][2] - a[0][2] * a[1][1])
    b[1][0] = t1 * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
    b[1][1] = t1 * (a[0][0] * a[2][2] - a[0][2] * a[2][0])
    b[1][2] = t1 * (a[0][2] * a[1][0] - a[0][0] * a[1][2])
    b[2][0] = t1 * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    b[2][1] = t1 * (a[0][1] * a[2][0] - a[0][0] * a[2][1])
    b[2][2] = t1 * (a[0][0] * a[1][1] - a[0][1] * a[1][0])
    
    return

#
# angular part of site-centered orbital
#  current implementation is for d-orbitals only
#
class Orbital:

  def __init__(self,coefs):
    self.pos=[]
    self.tri=[]
    self.coefs=coefs
    self.make()
  
  # real spherical harmonics 
  def Rlm(self,l,m,theta,phi):
    if l==2 and m==-2: return -math.sqrt(15.0/(16*math.pi))*math.sin(2*phi)*math.sin(theta)**2
    if l==2 and m==-1: return -math.sqrt(15.0/(16*math.pi))*math.sin(phi)*math.sin(2*theta)
    if l==2 and m==0: return math.sqrt(5/(64*math.pi))*(1+3*math.cos(2*theta))
    if l==2 and m==1: return -math.sqrt(15.0/(16*math.pi))*math.cos(phi)*math.sin(2*theta)
    if l==2 and m==2: return math.sqrt(15.0/(16*math.pi))*math.cos(2*phi)*math.sin(theta)**2
  
  def val(self,theta,phi):
    v=0
    for m in range(5): 
      v+=self.coefs[m]*self.Rlm(2,m-2,theta,phi)
    return v
    
  def make(self):
    raw_pos=[]
    raw_con=[]
    n=30
    for t in range(n):
      theta=math.pi*t/(n-1)
      for p in range(n):
        phi=2*math.pi*p/(n-1)
        v=5.5*self.val(theta,phi)
        x=v*math.sin(theta)*math.cos(phi)
        y=v*math.sin(theta)*math.sin(phi)
        z=v*math.cos(theta)
        raw_pos.append([x,y,z])
    for t in range(n):
      for p in range(n):
        i1=t
        i2=(t+1)%n
        j1=p
        j2=(p+1)%n
        n1=i1*n+j1
        n2=i1*n+j2
        n3=i2*n+j2
        n4=i2*n+j1
        raw_con.append([n1,n2,n3])
        raw_con.append([n1,n3,n4])
    # find equal positions
    eq_pos=[-1 for i in range(n*n)]
    l=0
    for i in range(n*n):
      if eq_pos[i]==-1:
        eq_pos[i]=l
        self.pos.append(raw_pos[i])
        for j in range(i+1,n*n):
          if abs(raw_pos[i][0]-raw_pos[j][0])<1e-10 and \
             abs(raw_pos[i][1]-raw_pos[j][1])<1e-10 and \
             abs(raw_pos[i][2]-raw_pos[j][2])<1e-10:
            eq_pos[j]=l
        l+=1
    npos=l
    # substitute positions in triangles by non-equal positions
    for i in range(2*n*n):
      raw_con[i][0]=eq_pos[raw_con[i][0]]
      raw_con[i][1]=eq_pos[raw_con[i][1]]
      raw_con[i][2]=eq_pos[raw_con[i][2]]

    eq_con=[-1 for i in range(2*n*n)]
    # mark degenerate triangles
    for i in range(2*n*n):
      if raw_con[i][0]==raw_con[i][1] or raw_con[i][0]==raw_con[i][2] or \
         raw_con[i][1]==raw_con[i][2]: eq_con[i]=-2
         
    # find equal triangles
    l=0
    for i in range(2*n*n):
      if eq_con[i]==-1:
        eq_con[i]=l
        self.tri.append(raw_con[i])
        for j in range(i+1,2*n*n):
          if raw_con[i][0]==raw_con[j][0] and raw_con[i][1]==raw_con[j][1] and \
             raw_con[i][2]==raw_con[j][2]: eq_con[j]=l
        l+=1     


#
# species-specific variables
#
class Species:
  
    def __init__(self, label):
        self.label = label
        self.R = 1.0
        self.color = [0.5, 0.5, 0.5]
        self.visible = True
    

# 
# atom-specific variables 
#
class Atom:

    def __init__(self, species, posc, posl):
        self.species = species
        self.posc = posc
        self.posl = posl
        self.nghbr = []
        self.orbital = 0
    
#
# geometry-specific variables
# 
class Geometry:

    def __init__(self):
        self.avec = []
        self.speciesList = {}
        self.atomList = []
        # read 'GEOMETRY.OUT'
        self.readGeometry()
        # make a list of nearest neighbours for each atom
        self.findNeighbours()
        # print basic info
        self.printGeometry()

    def readGeometry(self):
        
        fin = open("GEOMETRY.OUT","r")
        
        while True :
            
            line = fin.readline()
            if not line: break
            
            line = line.strip(" \n")
            
            if line == "avec":
                for i in range(3):
                    s1 = fin.readline().strip(" \n").split()
                    self.avec.append([float(s1[0]), float(s1[1]), float(s1[2])])
            
            if line == "atoms":
                # get number of species
                s1 = fin.readline().strip(" \n").split()
                nspecies = int(s1[0])
                
                # go over species
                for i in range(nspecies):
                    # construct label from species file name
                    s1 = fin.readline().strip(" \n").split()
                    label = s1[0][1:s1[0].find(".in")]
                    
                    # crate new species
                    sp = Species(label)
                    
                    # put species to the list
                    self.speciesList[label] = sp
                    
                    # get number of atoms for current species
                    s1 = fin.readline().strip(" \n").split()
                    natoms = int(s1[0])
                    # go over atoms
                    for j in range(natoms):
                        s1 = fin.readline().strip(" \n").split()
                        posl = [float(s1[0]), float(s1[1]), float(s1[2])]
                        posc = [0, 0, 0]
                        for l in range(3):
                            for x in range(3):
                                posc[x] += posl[l] * self.avec[l][x]
                    
                        # create new atom
                        self.atomList.append(Atom(sp, posc, posl))
        
        fin.close()
    
    def printGeometry(self):
        print "lattice vectors"
        print "  a1 : %12.6f %12.6f %12.6f"%(self.avec[0][0], self.avec[0][1], self.avec[0][2])
        print "  a2 : %12.6f %12.6f %12.6f"%(self.avec[1][0], self.avec[1][1], self.avec[1][2])
        print "  a3 : %12.6f %12.6f %12.6f"%(self.avec[2][0], self.avec[2][1], self.avec[2][2])
        print "atoms"
        for i in range(len(self.atomList)):
            print "%4i (%2s) at position %12.6f %12.6f %12.6f"%\
                (i, self.atomList[i].species.label, self.atomList[i].posc[0],\
                self.atomList[i].posc[1], self.atomList[i].posc[2])

    def findNeighbours(self):
        for iat in range(len(self.atomList)):
            xi = self.atomList[iat].posc
            nn = []
            # add nearest neigbours
            for jat in range(len(self.atomList)):
                xj = self.atomList[jat].posc
                for i1 in range(-4,5):
                    for i2 in range(-4,5):
                        for i3 in range(-4,5):
                            t = [0, 0, 0]
                            for x in range(3):
                                t[x] = i1 * self.avec[0][x] + i2 * self.avec[1][x] + i3 * self.avec[2][x]
                            r = [0, 0, 0]
                            for x in range(3):
                                r[x] = xj[x] + t[x] - xi[x]
                            d = math.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                            if (d <= 10.0):
                                nn.append([jat, r, d])
            # sort by distance
            for i in range(len(nn) - 1):
                for j in range(i+1, len(nn)):
                    if nn[j][2] < nn[i][2]:
                        nn[i], nn[j] = nn[j], nn[i]
            
            self.atomList[iat].nghbr = nn[:]


#
# cell (not necessarily primitive) with atoms and bonds
#
class Cell:
  
    def __init__(self, geometry, box):
        self.geometry = geometry
        self.box = box
        self.bonds = []
        self.atoms = []
        self.bondList = []
        return
      
    def hide(self, label):
        print " "
        print "hiding", label
        self.geometry.speciesList[label].visible = False
        return

    def atomSphere(self, label, color, R):
        self.geometry.speciesList[label].color = color
        self.geometry.speciesList[label].R = R
        return

    def bond(self, label1, label2, length, extend):
        self.bondList.append([label1, label2, length, extend])
        return

    def atomOrbital(self, ias, fname, iorb):
        fin = open(fname, "r")
        f1 = []
        for i in range(5):
            s1 = fin.readline().strip(" \n").split()
            if i == (iorb - 1):
                for j in range(5): 
                    f1.append(float(s1[j]))
        self.geometry.atomList[ias].orbital = Orbital(f1)
        return
      
    def write(self):
      self.fillBox()
      self.makeBonds()
      self.writeAtoms()
      self.writeBonds()
      #self.writeOrbitals()
      
    
#    def inBox(self, p, box):
#      n=[0,0,0]
#      a=box[0]
#      b=box[1]
#      c=box[2]
#      i=0
#      
#      n[0]=a[1]*b[2]-a[2]*b[1] 
#      n[1]=a[2]*b[0]-a[0]*b[2] 
#      n[2]=a[0]*b[1]-a[1]*b[0] 
#      d1=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]
#      d2=n[0]*(p[0]-c[0])+n[1]*(p[1]-c[1])+n[2]*(p[2]-c[2])
#      if cmp(d1,0)*cmp(d2,0)==-1 or abs(d1) < 1e-4 or abs(d2) < 1e-4: 
#        i+=1;
#
#      n[0]=a[1]*c[2]-a[2]*c[1] 
#      n[1]=a[2]*c[0]-a[0]*c[2] 
#      n[2]=a[0]*c[1]-a[1]*c[0] 
#      d1=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]
#      d2=n[0]*(p[0]-b[0])+n[1]*(p[1]-b[1])+n[2]*(p[2]-b[2])
#      if cmp(d1,0)*cmp(d2,0)==-1 or abs(d1) < 1e-4 or abs(d2) < 1e-4: 
#        i+=1;
#
#      n[0]=b[1]*c[2]-b[2]*c[1] 
#      n[1]=b[2]*c[0]-b[0]*c[2] 
#      n[2]=b[0]*c[1]-b[1]*c[0] 
#      d1=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]
#      d2=n[0]*(p[0]-a[0])+n[1]*(p[1]-a[1])+n[2]*(p[2]-a[2])
#      if cmp(d1,0)*cmp(d2,0)==-1 or abs(d1) < 1e-4 or abs(d2) < 1e-4: 
#        i+=1;
#        
#      if i==3: return True
#      else: return False

    def inBox(self, r, imv):
        # coorinates in the inits of box vectors
        rb = [0, 0, 0]
        for i in range(3):
            for j in range(3):
                rb[i] += imv[i][j] * r[j]

        if (rb[0] >= -0.5 and rb[0] <= 0.5) and \
           (rb[1] >= -0.5 and rb[1] <= 0.5) and \
           (rb[2] >= -0.5 and rb[2] <= 0.5): return True
        else: return False

    def fillBox(self):
        print " "
        print "populating the box"
        print "  box parameters"
        print "    center : %12.6f %12.6f %12.6f"%(self.box[0][0], self.box[0][1], self.box[0][2])
        print "        v1 : %12.6f %12.6f %12.6f"%(self.box[1][0], self.box[1][1], self.box[1][2])
        print "        v2 : %12.6f %12.6f %12.6f"%(self.box[2][0], self.box[2][1], self.box[2][2])
        print "        v3 : %12.6f %12.6f %12.6f"%(self.box[3][0], self.box[3][1], self.box[3][2])

        mv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        for i in range(3):
            for x in range(3):
                mv[x][i] = self.box[1+i][x]
          
        imv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        r3minv(mv, imv)
        
        for ias in range(len(self.geometry.atomList)):
            if self.geometry.atomList[ias].species.visible:
                for i1 in range(-4, 5):
                    for i2 in range(-4, 5):
                        for i3 in range(-4, 5):
                            # absolute position (position in the unit cell + translation)
                            r = [0, 0, 0]
                            for x in range(3):
                                r[x] = self.geometry.atomList[ias].posc[x] + \
                                       i1 * self.geometry.avec[0][x] + \
                                       i2 * self.geometry.avec[1][x] + \
                                       i3 * self.geometry.avec[2][x]
                            # position with respect to the center of the box
                            r0 = [0, 0, 0]
                            for x in range(3):
                                r0[x] = r[x] - self.box[0][x]

                            if self.inBox(r0, imv):
                                self.atoms.append([ias, r])
        return

    def writeAtoms(self):
        print " "
        print "writing ATOMS.dx"
        fout = open("ATOMS.dx", "w+")
        fout.write("object 1 class array type float rank 0 items %i data follows\n"%len(self.atoms))
        for i in range(len(self.atoms)):
            ias = self.atoms[i][0]
            fout.write("%f\n"%self.geometry.atomList[ias].species.R)
        fout.write("attribute \"dep\" string \"positions\"\n") 
        fout.write("#\n")
        fout.write("object 2 class array type float rank 1 shape 3 items %i data follows\n"%len(self.atoms))
        for i in range(len(self.atoms)):
            ias = self.atoms[i][0]
            color = self.geometry.atomList[ias].species.color
            fout.write("%f %f %f\n"%(color[0], color[1], color[2]))
        fout.write("attribute \"dep\" string \"positions\"\n") 
        fout.write("#\n")
        fout.write("object 3 class array type float rank 1 shape 3 items %i data follows\n"%len(self.atoms))
        for i in range(len(self.atoms)):
            ias = self.atoms[i][0]
            pos = self.atoms[i][1]
            fout.write("%f %f %f # %s\n"%(pos[0], pos[1], pos[2], self.geometry.atomList[ias].species.label))
        fout.write("attribute \"dep\" string \"positions\"\n") 
        fout.write("#\n")
        fout.write("object \"atoms\" class field\n")
        fout.write("component \"data\" value 1\n")
        fout.write("component \"colors\" value 2\n")
        fout.write("component \"positions\" value 3\n")
        fout.write("attribute \"name\" string \"cell\"")
        fout.close()
        return

    def index_in_atoms(self, r):
        for i in range(len(self.atoms)):
            if math.fabs(r[0] - self.atoms[i][1][0]) < 1e-10 and \
               math.fabs(r[1] - self.atoms[i][1][1]) < 1e-10 and \
               math.fabs(r[2] - self.atoms[i][1][2]) < 1e-10: return i
        return -1

    def makeBonds(self):
        for ibond in range(len(self.bondList)):
            lbl1 = self.bondList[ibond][0]
            lbl2 = self.bondList[ibond][1]
            length = self.bondList[ibond][2]
            extend = self.bondList[ibond][3]
            # go over all atoms in the box
            for i in range(len(self.atoms)):
                ias = self.atoms[i][0]
                if self.geometry.atomList[ias].species.label == lbl1:
                    # go over nearest neigbours of atom ias
                    for j in range(len(self.geometry.atomList[ias].nghbr)):
                        jas = self.geometry.atomList[ias].nghbr[j][0]
                        if (self.geometry.atomList[jas].species.label == lbl2) and \
                           (self.geometry.atomList[ias].nghbr[j][2] <= length):
                            # absolute position of neigbour: position of central atom + connecting vector
                            rj = [0, 0, 0]
                            for x in range(3):
                                rj[x] = self.atoms[i][1][x] + self.geometry.atomList[ias].nghbr[j][1][x]
                            # index of this neigbour in the list of atoms in the box
                            idx = self.index_in_atoms(rj)

                            if idx!=-1:
                                self.bonds.append([i, idx])
                            
                            elif extend:
                                self.atoms.append([jas, rj])
                                self.bonds.append([i, len(self.atoms)-1])
        return  

    def writeBonds(self):
        print " "
        print "writing BONDS.dx"
        fout = open("BONDS.dx","w+")
        fout.write("object 1 class array type float rank 1 shape 3 items %i data follows\n"%len(self.atoms))
        for i in range(len(self.atoms)):
            pos = self.atoms[i][1]
            fout.write("%f %f %f\n"%(pos[0], pos[1], pos[2]))
        fout.write("#\n")
        fout.write("object 2 class array type int rank 1 shape 2 items %i data follows\n"%len(self.bonds))
        for i in range(len(self.bonds)):
            fout.write("%i %i\n"%(self.bonds[i][0], self.bonds[i][1]))
        fout.write("attribute \"element type\" string \"lines\"\n")
        fout.write("attribute \"ref\" string \"positions\"\n")
        fout.write("#\n")
        fout.write("object \"atom_connect\" class field\n")
        fout.write("component \"positions\" value 1\n")
        fout.write("component \"connections\" value 2\n")
        fout.write("end\n")
        fout.close()
        return

    def writeOrbitals(self):
      print " "
      print "writing ORBITALS.dx"
      fout=open("ORBITALS.dx","w+")
      iorb=0
      for iat in range(len(self.atomList)):
        print self.atomList[iat].orbital
        if self.atomList[iat].orbital != 0:
          iorb+=1
          r0=self.atomList[iat].posc
          fout.write("object %i class array type float rank 1 shape 3 items %i data follows\n"%\
            ((iorb-1)*2+1,len(self.atomList[iat].orbital.pos)))
          for i in range(len(self.atomList[iat].orbital.pos)):
            r=[0,0,0]
            for x in range(3):
              r[x]=r0[x]+self.atomList[iat].orbital.pos[i][x]
            fout.write("%f %f %f\n"%(r[0],r[1],r[2]))
          fout.write("#\n")
          fout.write("object %i class array type int rank 1 shape 3 items %i data follows\n"%\
            ((iorb-1)*2+2,len(self.atomList[iat].orbital.tri)))  
          for i in range(len(self.atomList[iat].orbital.tri)):
            fout.write("%i %i %i\n"%(self.atomList[iat].orbital.tri[i][0],\
                                     self.atomList[iat].orbital.tri[i][1],\
                                     self.atomList[iat].orbital.tri[i][2]))
          fout.write("attribute \"ref\" string \"positions\"\n")
          fout.write("attribute \"element type\" string \"triangles\"\n")
          fout.write("attribute \"dep\" string \"connections\"\n")
          fout.write("#\n")  
          fout.write("object \"orbital%i\" class field\n"%iorb)
          fout.write("component \"positions\" value %i\n"%((iorb-1)*2+1))
          fout.write("component \"connections\" value %i\n"%((iorb-1)*2+2))
          fout.write("#\n")
      norb=iorb
      fout.write("object \"orbital\" class group\n")
      for iorb in range(norb):
        fout.write("  member %i value \"orbital%i\"\n"%(iorb,iorb+1))
      fout.write("end\n")
      fout.close()
      

#
# 
# 
print " "
print "GEOMTERY.OUT to OpenDX"
print " "
#
# get the geometry
#
geometry = Geometry()

#
# 3D box (center point + 3 non-collinear vectors) 
#  example: 
#   box=[[0,0,0],[10,0,0],[0,10,0],[0,0,10]] 
box = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]] #geometry.avec[0],geometry.avec[1],geometry.avec[2]]

#
# cell with user-defined shape
# 
cell = Cell(geometry, box)

# 
# cell.hide(label)
#  hides species with a given label
#  example:
#   cell.hide("Ca")
#cell.hide("Y")

#
# cell.atomSphere(label,color,radius)
#  sets [r,g,b] color and radius of species with a given label
#  example: red Mn sphere with radius 1.0
#   cell.atomSphere("Mn",[1.0,0.0,0.0],1.0)
cell.atomSphere("La", [0.0, 1.0, 0.0], 1.0)
cell.atomSphere("Cu", [1.0, 0.0, 0.0], 1.0)
cell.atomSphere("O", [0.0, 0.0, 1.0], 1.0)

#
# cell.bond(label1,label2,d,extend)
#  defines a bond with a maximum length 'd' from species 'label1' to species 'label2'
#  if extend is True, the bond can go outside the box
#  example: Mn-O bond
#   cell.bond("Mn","O",5,True)
cell.bond("Cu", "O", 5, True)

#
# cell.atomOrbital(j, file_name, i)
#  defines angular part of the site-centered orbital i for atom j; 
#  the orbital coefficients are taken from file file_name
#  example: 
#cell.atomOrbital(4, "Cu1_mtrx.txt", 1)

# 
# write to .dx files
#
cell.write()



