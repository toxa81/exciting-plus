#!/usr/bin/python

#
# linear response post processing 
#  
# Created June 2009 (AVK) 
#

import os
import sys
import getopt
import glob


class datFile:
  
  def __init__(self,name):
    self.name=name
    self.label=""
    self.title="title"
    self.qdir="dir"
    self.kmesh=[]
    self.emax=0.0
    self.de=0.0
    self.eta=0.0
    self.q=0.0
    self.ngv=0
    self.fxca=""
    self.fxctype=0
    self.data=[]
    self.parse(name)
    return
  
  def parse(self,name):
    lines=open(fileName,'r').readlines()
    for line in lines:
      if line.find(" : ")>-1:
        s2=line[line.index(":")+1:].strip()
        s1=s2.split()
      if line.find("k-mesh division")>-1:
        self.kmesh=[int(s1[0]),int(s1[1]),int(s1[2])]
      if line.find("maximum energy")>-1:
        self.emax=float(s1[0])
      if line.find("energy step")>-1:
        self.de=float(s1[0])
      if line.find("eta ")>-1:
        self.eta=float(s1[0])
      if line.find("q-vector length         [1/A]")>-1:
        self.q=float(s1[0])
      if line.find("G-vectors")>-1:
        self.ngv=int(s1[1])-int(s1[0])+1
      if line.find("label")>-1:
        self.label=s2
        print s2
      if line[0]!='#':
        t1=line.split()
        for i in range(len(t1)):
          t1[i]=float(t1[i])
        self.data.append(t1)
    return
  
  def plotFile(self,title,qdir):
    self.title=title.strip()
    self.qdir=qdir.strip()
    (name,extension)=os.path.splitext(self.name)
    name="plot__"+name
    plotname=name+".gnu"
    out=open(plotname,'w')
#    out.write("set terminal pdf dashed fsize 14 size 8.5,11\n")
    out.write("set terminal postscript portrait color\n")
    out.write("set output '| ps2pdf - "+name+".pdf'\n")
    q=("%6.3f"%self.q).strip()
    eta=("%6.3f"%self.eta).strip()
    str="set multiplot layout 2,2 title \"\\n"+self.title+", q="+q+" [1/A] || "+\
      self.qdir+"\\n ["+repr(self.kmesh[0])+"x"+repr(self.kmesh[1])+"x"+repr(self.kmesh[2])+"] k-mesh, "+\
      repr(self.ngv)+" G-vectors, eta="+eta+" eV"
    if self.fxca!="":
	  str=str+", fxcA="+self.fxca
    str=str+"\\n details : units are [eV], [1/eV/A^3], "+self.label+"\""
    out.write(str+"\n")
    out.write("set tmargin 2\n")
    out.write("set bmargin 0\n")
    out.write("set rmargin 1\n")
    out.write("set lmargin 1.5\n")
    out.write("LW=3\n")
### chiKS ###
#    out.write("set ylabel '[1/eV/nm^3]' offset 2,0\n")
#    out.write("unset xlabel\n")
#    out.write("set origin 0,0.44\n")
#    out.write("set size 0.46,0.44\n")
    str="plot [0:"+repr(self.emax)+"] "+"'"+self.name+"'"+\
      " using 1:($2)*1000 with lines title '-Re chiKS' lw LW, "+"'"+self.name+"'"+\
      " using 1:($3)*1000 with lines title '-Im chiKS' lw LW lt 1 lc rgb \""+"\x23"+"228b22\""
    out.write(str+"\n")
### chi ###
#    out.write("set ylabel '[1/eV/nm^3]' offset 2,0\n")
#    out.write("unset xlabel\n")
#    out.write("set origin 0.48,0.44\n")
#    out.write("set size 0.46,0.44\n")
    str="plot [0:"+repr(self.emax)+"] "+"'"+self.name+"'"+\
      " using 1:($5)*1000 with lines title '-Im chi' lw LW, "+"'"+self.name+"'"+\
      " using 1:($8)*1000 with lines title '-Im chi_scal' lw LW lc rgb \""+"\x23"+"FF0000\""
    out.write(str+"\n")
### epsilon_GqGq
#    out.write("set xlabel 'Energy [eV]'\n")
#    out.write("unset ylabel\n")
#    out.write("set origin 0,0\n")
#    out.write("set size 0.46,0.44\n") 
    str="plot [0:"+repr(self.emax)+"] "+"'"+self.name+"'"+\
      " using 1:11 with lines title '-Re eps_GqGq' lw LW, "+"'"+self.name+"'"+\
      " using 1:12 with lines title '-Im eps_GqGq' lw LW lt 1 lc rgb \""+"\x23"+"228b22\""
    out.write(str+"\n")
### epsilon_eff
#    out.write("set xlabel 'Energy [eV]'\n")
#    out.write("unset ylabel\n")
#    out.write("set origin 0.48,0\n")
#    out.write("set size 0.46,0.44\n")
    str="plot [0:"+repr(self.emax)+"] "+"'"+self.name+"'"+\
      " using 1:9 with lines title '-Re eps_eff' lw LW, "+"'"+self.name+"'"+\
      " using 1:10 with lines title '-Im eps_eff' lw LW lt 1 lc rgb \""+"\x23"+"228b22\""
    out.write("set label '"+self.name+"' font 'Times,12' at -74,9.2\n")
    out.write(str+"\n")
##   
    out.write("unset multiplot\n")

    out.close()
    os.system("gnuplot "+plotname)
    #os.system("rm "+plotname)
    


print " "
print "linear response post processing"
print " "

df=[]
for fileName in glob.glob("q*.dat"):
    df.append(datFile(fileName))
print "%i file(s) loaded\n"%len(df)

print("Input title : ")
title=sys.stdin.readline()
print("Input q-direction : ")
qdir=sys.stdin.readline()


for i in range(len(df)):
  df[i].plotFile(title,qdir)
