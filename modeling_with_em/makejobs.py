from sys import argv
import os
import shutil




#v                                                                                  xl penalty weight of gmm
jobs_table=[("template",  "modeling-1.1",  "modeling.py","n82.1.1","32"),
            ("template",  "modeling-1.2",  "modeling.py","n82.1.2","32")]



for n,job in enumerate(jobs_table):
  
  template_directory=job[0]
  directory_name=job[1]
  model_dir=os.getcwd()+"/"+directory_name
  script=job[2] #ex. bigj-ac-40s.py
  job_name=job[3]
  nprocesses=job[4]  
  output_file_name=model_dir
  error_file_name=model_dir
  
  info=''
  info+="# "+str(n)+" job name           "+job_name+"\n"
  info+="# "+str(n)+" template directory "+template_directory+"\n"
  info+="# "+str(n)+" final directory    "+directory_name+"\n"
  info+="# "+str(n)+" script name        "+script+"\n"
  info+="#########################################"+"\n"
  
  print info
  
  #shutil.copytree(template_directory, directory_name)
  os.mkdir(directory_name)
  shutil.copy(template_directory+"/"+script,directory_name+"/"+script)
  os.chdir(directory_name)
  
  jobscript='''
#!/bin/bash
#
#$ -S /bin/bash
#$ -l netapp=1G,scratch=1G
#$ -cwd
#$ -o ####OUTPUT####
#$ -e ####ERROR####
#$ -j y
#$ -l arch=linux-x64
#$ -l mem_free=2G
#$ -pe ompi ####NPROCESSES####
#$ -R yes
#$ -l netappsali=2G                  
#$ -l scrapp=2G
#$ -V
#$ -l h_rt=300:00:0.
#$ -t 1
#$ -N ####JOBNAME####

####INFO####
  
# load MPI modules
module load openmpi-1.6-nodlopen
module load /salilab/diva1/home/modules/sali-libraries
# IMP stuff

export IMP=/netapp/sali/pellarin/imp-050914/imp-fast-openmpi/setup_environment.sh

# write hostname and starting time 
hostname
date

# run the job
mpirun -np $NSLOTS $IMP python ####SCRIPT####

# done
date
'''

  jobscript=jobscript.replace("####OUTPUT####",output_file_name)
  jobscript=jobscript.replace("####ERROR####",error_file_name)
  jobscript=jobscript.replace("####JOBNAME####",job_name)
  jobscript=jobscript.replace("####SCRIPT####",script)
  jobscript=jobscript.replace("####NPROCESSES####",nprocesses)
  jobscript=jobscript.replace("####INFO####",info)
  jobscript=jobscript.replace("####MODELDIR####",directory_name)
  jobfile=open("job.sh","w")
  jobfile.write(jobscript)
  
  os.chdir("../")


