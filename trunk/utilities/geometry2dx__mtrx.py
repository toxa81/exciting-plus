#!/usr/bin/python

#
# GEOMTERY.OUT to OpenDX 
#  
# Created April 2009 (AVK) 
#

import math

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
  
  def __init__(self,label):
    self.label=label
    self.R=1.0
    self.color=[0.5,0.5,0.5]
    self.visible=True
    self.orbital=0
    

# 
# atom-specific variables 
#
class Atom:

  def __init__(self,species,posc,posl):
    self.species=species
    self.posc=posc
    self.posl=posl
    self.nghbr=[]
    
#
# geometry-specific variables
# 
class Geometry:

  def __init__(self):
    self.avec=[]
    self.speciesList={}
    self.atomList=[]
    # read 'GEOMETRY.OUT'
    self.readGeometry()
    # make a list of nearest neighbours for each atom
    self.findNeighbours()
    # print basic info
    self.printGeometry()

  def readGeometry(self):
    fin=open("GEOMETRY.OUT","r")
    while 1:
      line=fin.readline()
      if line=="": break
      line=line.strip(" \n")
      
      if line=="avec":
        for i in range(3):
          s1=fin.readline().strip(" \n").split()
          f1=[]
          for j in range(3):
            f1.append(float(s1[j]))
          self.avec.append(f1)
      
      if line=="atoms":
        # get number of species
        s1=fin.readline().strip(" \n").split()
        nspecies=int(s1[0])
        # go over species
        for i in range(nspecies):
          # construct label from species file name
          s1=fin.readline().strip(" \n").split()
          label=s1[0][1:s1[0].find(".in")]
          # crate new species
          sp=Species(label)
          # put species to the list
          self.speciesList[label]=sp
          # get number of atoms for current species
          s1=fin.readline().strip(" \n").split()
          natoms=int(s1[0])
          # go over atoms
          for j in range(natoms):
            s1=fin.readline().strip(" \n").split()
            posl=[]
            for x in range(3):
              posl.append(float(s1[x]))
            posc=[0,0,0]
            for l in range(3):
              for x in range(3):
                posc[x]+=posl[l]*self.avec[l][x]
            # create new atom
            self.atomList.append(Atom(sp,posc,posl))
    fin.close()
  
  def printGeometry(self):
    print "lattice vectors"
    print "  a1 : %12.6f %12.6f %12.6f"%(self.avec[0][0],self.avec[0][1],self.avec[0][2])
    print "  a2 : %12.6f %12.6f %12.6f"%(self.avec[1][0],self.avec[1][1],self.avec[1][2])
    print "  a3 : %12.6f %12.6f %12.6f"%(self.avec[2][0],self.avec[2][1],self.avec[2][2])
    print "atoms"
    for i in range(len(self.atomList)):
      print "%4i (%2s) at position %12.6f %12.6f %12.6f"%\
        (i,self.atomList[i].species.label,self.atomList[i].posc[0],\
        self.atomList[i].posc[1],self.atomList[i].posc[2])

  def findNeighbours(self):
    for iat in range(len(self.atomList)):
      xi=self.atomList[iat].posc
      nn=[]
      # add nearest neigbours
      for jat in range(len(self.atomList)):
        xj=self.atomList[jat].posc
        for i1 in range(-4,5):
          for i2 in range(-4,5):
            for i3 in range(-4,5):
              t=[0,0,0]
              for x in range(3):
                t[x]=i1*self.avec[0][x]+i2*self.avec[1][x]+\
                     i3*self.avec[2][x]
              r=[0,0,0]
              for x in range(3):
                r[x]=xj[x]+t[x]-xi[x]
              d=math.sqrt(r[0]**2+r[1]**2+r[2]**2)
              if (d<=10.0):
                nn.append([d,self.atomList[jat],r])
      # sort by distance
      for i in range(len(nn)-1):
        for j in range(i+1,len(nn)):
          if nn[j][0]<nn[i][0]:
            nn[i],nn[j]=nn[j],nn[i]
      
      self.atomList[iat].nghbr=nn[:]


#
# cell (not necessarily primitive) with atoms and bonds
#
class Cell:
  
  def __init__(self,geometry,box):
    self.geometry=geometry
    self.box=box
    self.bonds=[]
    self.atomList=[]
    self.bondList=[]
    
  def hide(self,label):
    print " "
    print "hiding",label
    self.geometry.speciesList[label].visible=False

  def atomSphere(self,label,color,R):
    self.geometry.speciesList[label].color=color
    self.geometry.speciesList[label].R=R

  def bond(self,label1,label2,length,extend):
    self.bonds.append([label1,label2,length,extend])
    return

  def atomOrbital(self,label,fname,row):
    fin=open(fname,"r")
    f1=[]
    for i in range(5):
      s1=fin.readline().strip(" \n").split()
      if i==(row-1):
        for j in range(5):
          f1.append(float(s1[j]))
    self.geometry.speciesList[label].orbital=Orbital(f1)
    
  def write(self):
    self.fillBox()
    self.makeBondList()
    self.writeAtoms()
    self.writeBonds()
    self.writeOrbitals()
    
  
  def inBox(self,p,box):
    n=[0,0,0]
    a=box[0]
    b=box[1]
    c=box[2]
    i=0
    
    n[0]=a[1]*b[2]-a[2]*b[1] 
    n[1]=a[2]*b[0]-a[0]*b[2] 
    n[2]=a[0]*b[1]-a[1]*b[0] 
    d1=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]
    d2=n[0]*(p[0]-c[0])+n[1]*(p[1]-c[1])+n[2]*(p[2]-c[2])
    if cmp(d1,0)*cmp(d2,0)==-1 or abs(d1) < 1e-4 or abs(d2) < 1e-4: 
      i+=1;

    n[0]=a[1]*c[2]-a[2]*c[1] 
    n[1]=a[2]*c[0]-a[0]*c[2] 
    n[2]=a[0]*c[1]-a[1]*c[0] 
    d1=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]
    d2=n[0]*(p[0]-b[0])+n[1]*(p[1]-b[1])+n[2]*(p[2]-b[2])
    if cmp(d1,0)*cmp(d2,0)==-1 or abs(d1) < 1e-4 or abs(d2) < 1e-4: 
      i+=1;

    n[0]=b[1]*c[2]-b[2]*c[1] 
    n[1]=b[2]*c[0]-b[0]*c[2] 
    n[2]=b[0]*c[1]-b[1]*c[0] 
    d1=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]
    d2=n[0]*(p[0]-a[0])+n[1]*(p[1]-a[1])+n[2]*(p[2]-a[2])
    if cmp(d1,0)*cmp(d2,0)==-1 or abs(d1) < 1e-4 or abs(d2) < 1e-4: 
      i+=1;
      
    if i==3: return True
    else: return False
  
  def fillBox(self):
    print " "
    print "populating the box"
    print "  box parameters"
    print "    center : %12.6f %12.6f %12.6f"%(self.box[0][0],self.box[0][1],self.box[0][2])
    print "        v1 : %12.6f %12.6f %12.6f"%(self.box[1][0],self.box[1][1],self.box[1][2])
    print "        v2 : %12.6f %12.6f %12.6f"%(self.box[2][0],self.box[2][1],self.box[2][2])
    print "        v3 : %12.6f %12.6f %12.6f"%(self.box[3][0],self.box[3][1],self.box[3][2])

    # find position of bottom-left-front corner of the box
    b0=[0,0,0]
    for x in range(3):
      b0[x]=self.box[0][x]-0.5*(self.box[1][x]+self.box[2][x]+self.box[3][x])
      
    for iat in range(len(self.geometry.atomList)):
      visible=self.geometry.atomList[iat].species.visible
      if visible:
        for i1 in range(-4,5):
          for i2 in range(-4,5):
            for i3 in range(-4,5):
              r=[0,0,0]
              r0=[0,0,0]
              t=[0,0,0]
              for x in range(3):
                t[x]=i1*self.geometry.avec[0][x]+i2*self.geometry.avec[1][x]+\
                  i3*self.geometry.avec[2][x]
              for x in range(3):
                r[x]=t[x]+self.geometry.atomList[iat].posc[x]
                r0[x]=r[x]-b0[x]
                
              if self.inBox(r0,[self.box[1],self.box[2],self.box[3]]):
                tmp=Atom(self.geometry.atomList[iat].species,r,[0,0,0])
                tmp.nghbr=self.geometry.atomList[iat].nghbr
                self.atomList.append(tmp)

  def writeAtoms(self):
    print " "
    print "writing ATOMS.dx"
    fout=open("ATOMS.dx","w+")
    fout.write("object 1 class array type float rank 0 items %i data follows\n"%len(self.atomList))
    for i in range(len(self.atomList)):
      fout.write("%f\n"%self.atomList[i].species.R)
    fout.write("attribute \"dep\" string \"positions\"\n") 
    fout.write("#\n")
    fout.write("object 2 class array type float rank 1 shape 3 items %i data follows\n"%len(self.atomList))
    for i in range(len(self.atomList)):
      color=self.atomList[i].species.color
      fout.write("%f %f %f\n"%(color[0],color[1],color[2]))
    fout.write("attribute \"dep\" string \"positions\"\n") 
    fout.write("#\n")
    fout.write("object 3 class array type float rank 1 shape 3 items %i data follows\n"%len(self.atomList))
    for i in range(len(self.atomList)):
      pos=self.atomList[i].posc
      fout.write("%f %f %f # %s\n"%(pos[0],pos[1],pos[2],self.atomList[i].species.label))
    fout.write("attribute \"dep\" string \"positions\"\n") 
    fout.write("#\n")
    fout.write("object \"atoms\" class field\n")
    fout.write("component \"data\" value 1\n")
    fout.write("component \"colors\" value 2\n")
    fout.write("component \"positions\" value 3\n")
    fout.write("attribute \"name\" string \"cell\"")
    fout.close() 

  def inAtomList(self,r):
    for i in range(len(self.atomList)):
      if abs(r[0]-self.atomList[i].posc[0])<1e-10 and \
         abs(r[1]-self.atomList[i].posc[1])<1e-10 and \
         abs(r[2]-self.atomList[i].posc[2])<1e-10: return i
    return -1

  def makeBondList(self):
    for ibond in range(len(self.bonds)):
      lbl1=self.bonds[ibond][0]
      lbl2=self.bonds[ibond][1]
      length=self.bonds[ibond][2]
      extend=self.bonds[ibond][3]
      n0=len(self.atomList)
      # go over all atoms in the box
      for iat in range(n0):
        if self.atomList[iat].species.label==lbl1:
          # go over nearest neigbours of atom iat
          for jat in range(len(self.atomList[iat].nghbr)):
            if self.atomList[iat].nghbr[jat][1].species.label==lbl2 and\
               self.atomList[iat].nghbr[jat][0]<=length:
              r=[0,0,0]
              for x in range(3):
                r[x]=self.atomList[iat].posc[x]+self.atomList[iat].nghbr[jat][2][x]
              idx=self.inAtomList(r)
              if idx!=-1:
                self.bondList.append([iat,idx])
              elif extend:
                tmp=Atom(self.atomList[iat].nghbr[jat][1].species,r,[0,0,0])
                self.atomList.append(tmp)
                self.bondList.append([iat,len(self.atomList)-1])
    return  

  def writeBonds(self):
    print " "
    print "writing BONDS.dx"
    fout=open("BONDS.dx","w+")
    fout.write("object 1 class array type float rank 1 shape 3 items %i data follows\n"%len(self.atomList))
    for i in range(len(self.atomList)):
      pos=self.atomList[i].posc
      fout.write("%f %f %f\n"%(pos[0],pos[1],pos[2]))
    fout.write("#\n")
    fout.write("object 2 class array type int rank 1 shape 2 items %i data follows\n"%len(self.bondList))
    for i in range(len(self.bondList)):
      fout.write("%i %i\n"%(self.bondList[i][0],self.bondList[i][1]))
    fout.write("attribute \"element type\" string \"lines\"\n")
    fout.write("attribute \"ref\" string \"positions\"\n")
    fout.write("#\n")
    fout.write("object \"atom_connect\" class field\n")
    fout.write("component \"positions\" value 1\n")
    fout.write("component \"connections\" value 2\n")
    fout.write("end\n")
    fout.close()

  def writeOrbitals(self):
    print " "
    print "writing ORBITALS.dx"
    fout=open("ORBITALS.dx","w+")
    iorb=0
    for iat in range(len(self.atomList)):
      if self.atomList[iat].species.orbital!=0:
        iorb+=1
        r0=self.atomList[iat].posc
        fout.write("object %i class array type float rank 1 shape 3 items %i data follows\n"%\
          ((iorb-1)*2+1,len(self.atomList[iat].species.orbital.pos)))
        for i in range(len(self.atomList[iat].species.orbital.pos)):
          r=[0,0,0]
          for x in range(3):
            r[x]=r0[x]+self.atomList[iat].species.orbital.pos[i][x]
          fout.write("%f %f %f\n"%(r[0],r[1],r[2]))
        fout.write("#\n")
        fout.write("object %i class array type int rank 1 shape 3 items %i data follows\n"%\
          ((iorb-1)*2+2,len(self.atomList[iat].species.orbital.tri)))  
        for i in range(len(self.atomList[iat].species.orbital.tri)):
          fout.write("%i %i %i\n"%(self.atomList[iat].species.orbital.tri[i][0],\
                                   self.atomList[iat].species.orbital.tri[i][1],\
                                   self.atomList[iat].species.orbital.tri[i][2]))
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
geometry=Geometry()

#
# 3D box (center point + 3 non-collinear vectors) 
#  example: 
#   box=[[0,0,0],[10,0,0],[0,10,0],[0,0,10]] 
box=[[0,0,0],[2,0,0],[0,2,0],[0,0,-30]] #geometry.avec[0],geometry.avec[1],geometry.avec[2]]

#
# cell with user-defined shape
# 
cell=Cell(geometry,box)

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
cell.atomSphere("V1a",[1.0,0.0,0.0],1.0)
cell.atomSphere("V1b",[1.0,0.0,0.0],1.0)
cell.atomSphere("V2",[1.0,0.0,0.0],1.0)
cell.atomSphere("O",[0.0,0.0,1.0],0.7)

#
# cell.bond(label1,label2,d,extend)
#  defines a bond with a maximum length 'd' from species 'label1' to species 'label2'
#  if extend is True, the bond can go outside the box
#  example: Mn-O bond
#   cell.bond("Mn","O",5,True)
#cell.bond("V","O",5,True)
cell.bond("V1a","V1b",5.24,True)

#
# cell.atomOrbital(label,coefs)
#  defines angular part of the site-centered orbital for a species 'label'
#  example: 3z^2-r^2 orbital
#    cell.atomOrbital("Mn",[0,0,1,0,0])
#cell.atomOrbital("V1","V1.txt",3)

# 
# write to .dx files
#
cell.write()



