import os,sys,argparse
import fileinput,time
import glob
import re,bisect
import pandas as pd
import numpy as np
from operator import itemgetter
import read_write_file 
import association_with_regions


def read_summits_from_narrowpeak_bed(infile,species,expand,write_out=None,end_col=2,region_sorted=True):
    chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
         'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
         'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
    regions = {}
    with open(infile,'r') as inf:
        line = inf.readline()
        while line:
            #print(line)
            if not re.match("#",line):
                #print(line)
                sline = line.strip().split()
                chrom = sline[0]
                start = int(sline[1]) + int(sline[-1])
                end = start#int(sline[end_col])
                fold_enrichment = float(sline[6])
                if chrom in chroms and fold_enrichment>=4:
                    if chrom not in regions:
                        regions[chrom]=[]
                    regions[chrom].append([max(0,start-expand),end+expand])
                else:
                    pass                   
            line = inf.readline()
    if region_sorted:
        for chrom in regions:
            regions[chrom] = sorted(regions[chrom],key = itemgetter(0),reverse=False)
    if write_out:
        with open(write_out,'w') as outf:
            for chrom in regions:
                for region in regions[chrom]:
                    outf.write('{}\t{}\t{}\n'.format(chrom,region[0],region[1]))  
    return regions

def get_overlap_info(islands,regions):
    # for each island in islands[chrom], check if 1-overlap 0-non-overlap with regions[chrom]
    island_dict = {}
    for chrom in islands:
        if chrom not in regions:
            for island in islands[chrom]:
                #print(island);exit(0)
                ID = island[-1].split('\t')[0]
                island_dict[ID] = 0
        else:
            # check the regions in regions[chrom] are non-overlap and sorted
            # as use the summit info, all summits should not overlapped, expand < 75
            regionlist = []
            for region in regions[chrom]:
                regionlist.extend(region)
            for num in range(len(regionlist)-1):
                assert regionlist[num]<regionlist[num+1]
                            
            # for each island, check if overlapped with regions[chrom]
            for island in islands[chrom]:
                ID = island[-1].split('\t')[0]
                assert island[0]<=island[1]
                s = bisect.bisect_left(regionlist,island[0])
                e = bisect.bisect_right(regionlist,island[1])
                if s==e and s%2==0:
                    island_dict[ID] = 0
                else:
                    island_dict[ID] = 1
    return island_dict


def get_overlap_info_matrix(df,islands,name,species):
    
    #df = df.fillna('NA')  
    overlap_df = pd.DataFrame()
    for data_index in df.index:
        peak_file = '{}_peaks.narrowPeak'.format(data_index)
        region1 = read_summits_from_narrowpeak_bed(peak_file,species,expand=25,end_col=1)
        overlap_df[data_index]=pd.Series(get_overlap_info(islands,region1))
#         print(overlap_df);exit(0)
    overlap_df.to_csv(name,sep='\t')

    

def main():   

    CTCF_collection_file='samples.csv'
    CTCF_collection_df = pd.read_csv(CTCF_collection_file,index_col=0)
    
    species='hg38'
    union_CTCF='CTCF_union_summits_fe4_width_150_sorted.bed'
    all_regions = association_with_regions.read_regions_from_bed_not_sep_strand(union_CTCF,species)
    get_overlap_info_matrix(CTCF_collection_df,all_regions,'union_summits_EachDataOverlapInfo_150bp.csv',species)


            
            
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of bed fromat, union all the overlapping regions', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)

    main()
