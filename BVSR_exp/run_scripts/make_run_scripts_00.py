#!/usr/bin/env python

'''make a run script for each subset and output a qsub file'''
import string
qsubfile = open('../qsub.txt','w')
prescript = '00_bvsr'

for pop in ['AFA','CAU','HIS']:
#    for i in range(1,23):
    for i in range(1,5):
        newi = str(i)
        outfilename = 'run_' + prescript + '_' + pop + '_' + newi + '.sh'
        outfile = open(outfilename,'w')
        output = '''#!/bin/bash
#PBS -N R.bvsr.''' + pop + newi + '''\n#PBS -S /bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

'''
        outfile.write(output)
        outfile.write('time R --no-save < ' + prescript + '.r --args ' + newi + ' ' + pop + '\n')

        qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')
