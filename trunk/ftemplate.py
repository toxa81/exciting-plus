# ftemplate - Fortran Templates  

import sys

tidx=0

for line in sys.stdin.readlines():
  sline=line.strip()
# default: write the line
  wline=1
# scan for python code outside of a template
  if sline.find("@python")==0 and tidx==0:
    wline=0
    py_code=sline.replace("@python","").strip()
# execute the code
    exec(py_code)
# we are at the beginning of a template
  if sline.find("@template")==0 and sline.find("begin")!=-1: 
    wline=0
    tidx=1
    template=[]
# we are in the template
  if tidx==1:
    wline=0
    template.append(line)
# now we have the full template; let's unroll it
  if sline.find("@template")==0 and sline.find("end")!=-1: 
    wline=0
    tidx=0 
# remove first and last lines which are '@template begin' and '@template end'    
    template.pop()
    template.pop(0)
# scan for variables, body and initial python code
    variables=[]
    body=[]
    for l in template:
      sl=l.strip()
      if sl.find("@template")==0 and sl.find("variable")!=1:
        varname=sl.replace("@template","").replace("variable","").strip()
        variables.append(varname)
      if sl.find("@")!=0: 
        body.append(l)
      if sl.find("@python")==0:
        py_code_initial=sl.replace("@python","").strip()
# main trick: compose python code from initial code and templated text 
    tstr=""
    for l in body:
      tstr=tstr+l
# replace templates with variables    
    s="tstr"
    for v in variables:
      s=s+".replace(\"#"+v+"\","+v+")"
    exec(py_code_initial+" sys.stdout.write("+s+")");
# write normal line      
  if wline==1:
    sys.stdout.write(line)

