#!/usr/bin/python

import os
import glob

title="V2O3 (PM phase)"
qdir="[001]"


def plotFile(fileName):
    lines=open(fileName,'r').readlines()
    for line in lines:
        if line.find("k-mesh division")>-1:
            t1=line.split()
            kx=t1[4]
            ky=t1[5]
            kz=t1[6]
            print "kx,ky,kz : ", kx,ky,kz
        if line.find("q-vector length         [1/A]")>-1:
            t1=line.split()
            q="%10.3f"%float(t1[5])
            q=q.strip()
            print "q : ", q
        if line.find("G-vectors")>-1:
            t1=line.split()
            ngv=t1[4]
            print "ngv : ",ngv
        if line.find("eta            [eV]")>-1:
            t1=line.split()
            eta="%10.2f"%float(t1[4])
            eta=eta.strip()
            print "eta : ",eta
        if line.find("fxc A")>-1:
            t1=line.split()
            fxcA="%10.2f"%float(t1[4])
            fxcA=fxcA.strip()
            print "fxc A : ",fxcA
            
    (name,extension)=os.path.splitext(fileName)
    name="plot__"+name
    plotname=name+".gnu"
    out=open(plotname,'w')
    out.write("set terminal pdf dashed fsize 14 size 8.5,11\n")
    out.write("set output '"+name+".pdf'\n")
    str="set multiplot layout 2,2 title \"\\n"+title+", q="+q+" [1/A] || "+\
      qdir+"\\n ["+kx+"x"+ky+"x"+kz+"] k-mesh, "+ngv+" G-vectors, eta="+eta+" eV\"" 
    out.write(str+"\n")
    out.write("set tmargin 0\n")
    out.write("set bmargin 4\n")
    out.write("set rmargin 0\n")
    out.write("set lmargin 10\n")
    out.write("LW=6\n")
### chiKS ###
    out.write("set ylabel '[1/eV/nm^3]' offset 2,0\n")
    out.write("unset xlabel\n")
    out.write("set origin 0,0.44\n")
    out.write("set size 0.46,0.44\n")
    str="plot [0:15] "+"'"+fileName+"'"+\
      " using 1:($2)*1000 with lines title '-Re chiKS' lw LW, "+"'"+fileName+"'"+\
      " using 1:($3)*1000 with lines title '-Im chiKS' lw LW lt 1 lc rgb \""+"\x23"+"228b22\""
    out.write(str+"\n")
### chi ###
    out.write("set ylabel '[1/eV/nm^3]' offset 2,0\n")
    out.write("unset xlabel\n")
    out.write("set origin 0.48,0.44\n")
    out.write("set size 0.46,0.44\n")
    str="plot [0:15] "+"'"+fileName+"'"+\
      " using 1:($5)*1000 with lines title '-Im chi' lw LW, "+"'"+fileName+"'"+\
      " using 1:($8)*1000 with lines title '-Im chi_scal' lw LW lc rgb \""+"\x23"+"FF0000\""
    out.write(str+"\n")
### epsilon_GqGq
    out.write("set xlabel 'Energy [eV]'\n")
    out.write("unset ylabel\n")
    out.write("set origin 0,0\n")
    out.write("set size 0.46,0.44\n") 
    str="plot [0:15] "+"'"+fileName+"'"+\
      " using 1:11 with lines title '-Re eps_GqGq' lw LW, "+"'"+fileName+"'"+\
      " using 1:12 with lines title '-Im eps_GqGq' lw LW lt 1 lc rgb \""+"\x23"+"228b22\""
    out.write(str+"\n")
### epsilon_eff
    out.write("set xlabel 'Energy [eV]'\n")
    out.write("unset ylabel\n")
    out.write("set origin 0.48,0\n")
    out.write("set size 0.46,0.44\n")
    str="plot [0:15] "+"'"+fileName+"'"+\
      " using 1:9 with lines title '-Re eps_eff' lw LW, "+"'"+fileName+"'"+\
      " using 1:10 with lines title '-Im eps_eff' lw LW lt 1 lc rgb \""+"\x23"+"228b22\""
    out.write(str+"\n")
##   
    out.write("unset multiplot\n")
    out.close()
    os.system("gnuplot "+plotname)
    os.system("rm "+plotname)
for fileName in glob.glob("q*.dat"):
    plotFile(fileName)
