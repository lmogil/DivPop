#!/usr/bin/env python

'''make a run script for each subset and output a qsub file'''
import string
qsubfile = open('../qsub.txt','w')
prescript = '08_est_gen_corr_simu-exp.r'
for i in range(1,10):
    newi = str(i)
    outfilename = 'run_' + prescript + '_' + newi + '.sh'
    outfile = open(outfilename,'w')
    output = '''#!/bin/bash
#PBS -N gcta.simu-exp''' + newi + '''\n#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
PATH=$PBS_O_PATH
cd $PBS_O_WORKDIR

'''
    outfile.write(output)
    for pop in ['AFA CAU','CAU HIS','AFA HIS']:
        outfile.write('time R --no-save < ' + prescript + ' --args ' + pop + ' ' + newi + '\n')

    qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')
