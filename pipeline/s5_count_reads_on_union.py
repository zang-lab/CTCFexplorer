import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from operator import itemgetter


def main():               
    
    sample_file='samples.csv'
    gsm_collection = pd.read_csv(sample_file)['GSM'].tolist()
    union_CTCF='CTCF_union_summits_fe4_width_150_sorted.bed'
    
    slurm_dir='run_slurm'
    os.makedirs(slurm_dir,exist_ok=True)
    RPKM_dir='RPKM_csv'
    os.makedirs(RPKM_dir,exist_ok=True)
    
    for gsmID in gsm_collection:
        bam_file = '{}_clean.bam'.format(gsmID)
        # check existence
        if os.path.isfile('{}/{}.csv'.format(RPKM_dir,gsmID)):
            pass
        else:
             with open('{}/{}.slurm'.format(slurm_dir,gsmID),'w') as slurmout:
                slurmout.write('''#!/bin/bash\n
#SBATCH -n 1
#SBATCH --mem=100G
#SBATCH -t 6:00:00
#SBATCH -p standard 
''')
                #######################
                # REMEMBER TO CHANGE THE MEMORY !!!
                #########################
            
                slurmout.write('#SBATCH -o {}/{}.out\n\n'.format(slurm_dir,gsmID))
                slurmout.write('python utils.count_reads_on_regions.py -i {} -b {} -f bam -s hg38 -o {}/{}.csv'.format(union_CTCF,bam_file,RPKM_dir,gsmID))

    outdir = 'signal_on_union'
    indir = os.path.join(outdir, 'RPKM_csv')

    # Collect all CSV files
    infiles = glob.glob(os.path.join(indir, '*.csv'))

    # Read and concatenate data
    df_all = pd.concat([pd.read_csv(f, sep='\t', header=None, index_col=0, names=[os.path.splitext(os.path.basename(f))[0]]) for f in infiles], axis=1)

    # Round values and save
    df_all.round(2).to_csv(os.path.join(outdir, 'union_RPKM.csv'), sep='\t')


         
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-a', '--infile1', action = 'store', type = str,dest = 'infile1', help = 'input file to be compared/separated', metavar = '<file>')
    #parser.add_argument('-b', '--infile2', action = 'store', type = str,dest = 'infile2', help = 'input file to be compared as basic', metavar = '<file>')
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main()
