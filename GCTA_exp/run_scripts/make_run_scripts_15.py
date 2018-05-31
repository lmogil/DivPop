#!/usr/bin/env python

'''make a run script for each subset and output a qsub file'''
import string
qsubfile = open('../qsub.txt','w')
prescript = '15_est_gen_corr_simu-exp_lauren_imp.r'
for i in range(10):
    for pop in ['AFA CAU','CAU HIS','AFA HIS']:
        newi = str(i)
        outfilename = 'run_' + prescript + '_' + pop[0] + pop[4] + newi + '.sh'
        outfile = open(outfilename,'w')
        output = '''#!/bin/bash
#PBS -N gcta.simu-exp''' + newi + '''\n#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
PATH=$PBS_O_PATH
cd $PBS_O_WORKDIR

'''
        outfile.write(output)
        outfile.write('time R --no-save < ' + prescript + ' --args ' + pop + ' ' + newi + '\n')

        qsubfile.write('qsub run_scripts/' + outfilename + '\nsleep 3\n')
