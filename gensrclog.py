import commands

svninfo=commands.getoutput("svn info") 
svndiff=commands.getoutput("svn diff")

fout=open("srclog.f90","w")
fout.write("!Warning: this is automatically generated file\n")
fout.write("subroutine srclog\n")
fout.write("implicit none\n")
fout.write("open(190,file=\"source-code.log\",status=\"replace\",form=\"formatted\")\n")
fout.write("write(190,'(\"Output of ''svn info''\")')\n")
t=svninfo.splitlines()
for s in t:
  fout.write("write(190,'(\""+s+"\")')\n")
fout.write("write(190,*)\n")  
fout.write("write(190,'(\"Output of ''svn diff''\")')\n")
t=svndiff.splitlines()
for s in t:
  s=s.replace("'","''")
  s=s.replace("\"","\"\"")
  fout.write("write(190,'(\""+s+"\")')\n")
fout.write("close(190)\n")  
fout.write("return\n")
fout.write("end\n")

