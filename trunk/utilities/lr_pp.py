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
      if line.find("fxc A")>-1:
        self.fxca=s1[0]
      if line[0]!='#':
        t1=line.split()
        for i in range(len(t1)):
          t1[i]=float(t1[i])
        self.data.append(t1)
    return
  
  def plotFile(self,title,qdir):
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
    str=str+"\\n details : units are [eV], [1/eV/nm^3], "+self.label+"\""
    out.write(str+"\n")
    out.write("set tmargin 2\n")
    out.write("set bmargin 0\n")
    out.write("set rmargin 1\n")
    out.write("set lmargin 1.5\n")
    out.write("LW=2\n")
#    self.emax=30
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
      " using 1:($8)*1000 with lines title '-Im chi_scal' lw LW lc rgb \""+"\x23"+"0000FF\""
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
#    out.write("set label '"+self.name+"' font 'Times,12' at -74,9.2\n")
    out.write(str+"\n")
##   
    out.write("unset multiplot\n")

    out.close()
    os.system("gnuplot "+plotname)
    os.system("rm "+plotname)
    
def plotFile2(title,qdir,df1,df2):
  print "Input label1 for file %s : "%df1.name  
  label1=sys.stdin.readline().strip()
  print "Input label2 for file %s : "%df2.name
  label2=sys.stdin.readline().strip()
  (name1,extension1)=os.path.splitext(df1.name)
  (name2,extension2)=os.path.splitext(df2.name)
  name="plot__"+name1+"_vs_"+name2
  plotname=name+".gnu"
  out=open(plotname,'w')
  out.write("set terminal postscript portrait color \"Helvetica\" 11\n")
  out.write("set output '| ps2pdf - "+name+".pdf'\n")
# make a title
  line1=title+", q-dir : "+qdir
  q=("%6.3f"%df1.q).strip()
  line2=label1+" : q="+q+" [1/A]"
  if df1.fxca!="":
    line2=line2+", fxcA="+df1.fxca
  q=("%6.3f"%df2.q).strip()
  line3=label2+" : q="+q+" [1/A]"
  if df2.fxca!="":
    line3=line3+", fxcA="+df2.fxca
  str="set multiplot layout 2,2 title \""+line1+"\\n"+line2+"\\n"+line3
  str=str+"\\n details : units are [eV], [1/eV/nm^3]\""
  out.write(str+"\n")
  out.write("set tmargin 2\n")
  out.write("set bmargin 0\n")
  out.write("set rmargin 1\n")
  out.write("set lmargin 1.5\n")
  out.write("LW=2\n")
  print "Input maximum energy : "
  emax=float(sys.stdin.readline().strip())
# nice green color is 228b22
### chiKS ###
#    out.write("set origin 0,0.44\n")
#    out.write("set size 0.46,0.44\n")
  str="plot [0:"+repr(emax)+"] "+\
    "'"+df1.name+"'"+\
    " using 1:($2)*1000 with lines title '-Re chiKS("+label1+")' lw LW lt 2 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df1.name+"'"+\
    " using 1:($3)*1000 with lines title '-Im chiKS("+label1+")' lw LW lt 1 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df2.name+"'"+\
    " using 1:($2)*1000 with lines title '-Re chiKS("+label2+")' lw LW lt 2 lc rgb \""+"\x23"+"0000FF\","+\
    "'"+df2.name+"'"+\
    " using 1:($3)*1000 with lines title '-Im chiKS("+label2+")' lw LW lt 1 lc rgb \""+"\x23"+"0000FF\""
  out.write(str+"\n")
### chi ###
#    out.write("set origin 0.48,0.44\n")
#    out.write("set size 0.46,0.44\n")
  str="plot [0:"+repr(emax)+"] "+\
    "'"+df1.name+"'"+\
    " using 1:($5)*1000 with lines title '-Im chi("+label1+")' lw LW lt 1 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df1.name+"'"+\
    " using 1:($8)*1000 with lines title '-Im chi_scal("+label1+")' lw LW lt 2 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df2.name+"'"+\
    " using 1:($5)*1000 with lines title '-Im chi("+label2+")' lw LW lt 1 lc rgb \""+"\x23"+"0000FF\","+\
    "'"+df2.name+"'"+\
    " using 1:($8)*1000 with lines title '-Im chi_scal("+label2+")' lw LW lt 2 lc rgb \""+"\x23"+"0000FF\""
  out.write(str+"\n")
### epsilon_GqGq
#    out.write("set origin 0,0\n")
#    out.write("set size 0.46,0.44\n") 
  str="plot [0:"+repr(emax)+"] "+\
    "'"+df1.name+"'"+\
    " using 1:11 with lines title '-Re eps_GqGq("+label1+")' lw LW lt 2 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df1.name+"'"+\
    " using 1:12 with lines title '-Im eps_GqGq("+label1+")' lw LW lt 1 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df2.name+"'"+\
    " using 1:11 with lines title '-Re eps_GqGq("+label2+")' lw LW lt 2 lc rgb \""+"\x23"+"0000FF\","+\
    "'"+df2.name+"'"+\
    " using 1:12 with lines title '-Im eps_GqGq("+label2+")' lw LW lt 1 lc rgb \""+"\x23"+"0000FF\""
  out.write(str+"\n")
### epsilon_eff
#    out.write("set origin 0.48,0\n")
#    out.write("set size 0.46,0.44\n")
  str="plot [0:"+repr(emax)+"] "+\
    "'"+df1.name+"'"+\
    " using 1:9 with lines title '-Re eps_eff("+label1+")' lw LW lt 2 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df1.name+"'"+\
    " using 1:10 with lines title '-Im eps_eff("+label1+")' lw LW lt 1 lc rgb \""+"\x23"+"FF0000\","+\
    "'"+df2.name+"'"+\
    " using 1:9 with lines title '-Re eps_eff("+label2+")' lw LW lt 2 lc rgb \""+"\x23"+"0000FF\","+\
    "'"+df2.name+"'"+\
    " using 1:10 with lines title '-Im eps_eff("+label2+")' lw LW lt 1 lc rgb \""+"\x23"+"0000FF\""
  out.write(str+"\n")
  out.write("unset multiplot\n")

  out.close()
  os.system("gnuplot "+plotname)
  os.system("rm "+plotname)

  return

def plotFile3(df):
  print "Input column to plot (1 is energy) : "  
  icol=int(sys.stdin.readline().strip())-1
  print "Input scale factor : "
  scale=float(sys.stdin.readline().strip())
  print "Scan through : "
  print "  1 : q-vectors"
  print "  2 : fxcA parameter"
  iscan=int(sys.stdin.readline().strip())
  
  
  out=open("plot.dat","w")
  
  for j in range(len(df[0].data)):    
    t1=[]
    for i in range(len(df)):
      if iscan==1:
        a1=float(df[i].q)
      if iscan==2:
        a1=float(df[i].fxca)
      a2=df[i].data[j][icol]*scale
      t1.append([a1,a2])
    
    for i1 in range(len(df)-1):
      for i2 in range(i1+1,len(df)):
        if t1[i2][0] < t1[i1][0]:
          t1[i1],t1[i2]=t1[i2],t1[i1]
    
    for i in range(len(df)):
      out.write("%f %f %f\n"%(df[0].data[j][0],t1[i][0],t1[i][1]))
    
    out.write(" \n")
  
  out.close()
  
  out=open("plot.gnu","w")
  out.write("set terminal postscript landscape color \"Helvetica\" 11\n")
#  out.write("set output '| ps2pdf - plot.pdf'\n")
  out.write("set output 'plot.ps'\n")
  out.write("set pm3d map \n")
#  out.write("set palette rgbformulae 22,13,-31\n")
# grayscale
#  out.write("set palette defined ( 0 1 1 1, 1 0 0 0 )\n")
#  out.write("set view 30,300\n")
  out.write("splot 'plot.dat' with pm3d\n")
  out.close()
  os.system("gnuplot plot.gnu")
  os.system("rm plot.gnu")
  os.system("rm plot.dat")
  



print " "
print "linear response post processing"
print " "

df=[]
for fileName in glob.glob("q*.dat"):
    df.append(datFile(fileName))
print "%i file(s) loaded\n"%len(df)
print "Index of files : "
for i in range(len(df)):
  print "%i  :  %s"%(i,df[i].name)
print " "

print "Available plotting modes : "
print "  1 : batch plot of all loaded files"
print "  2 : comparison plot of two files"
print "  3 : scan-plot of all loaded filed"  
print "Input your choice : "
str=sys.stdin.readline().strip()
if not (str=="1" or str=="2" or str=="3"):
  print "Wrong input!"
  exit()
 
print("Input title : ")
title=sys.stdin.readline().strip()
print("Input q-direction : ")
qdir=sys.stdin.readline().strip()

if str=="1":
  for i in range(len(df)):
    df[i].plotFile(title,qdir)

if str=="2":
  print "Input index of file1 : "  
  i1=int(sys.stdin.readline().strip())
  print "Input index of file2 : "  
  i2=int(sys.stdin.readline().strip())
  plotFile2(title,qdir,df[i1],df[i2])

if str=="3":
  plotFile3(df)
  