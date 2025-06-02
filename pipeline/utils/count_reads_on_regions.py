# Time-stamp: <2017-08-10>
'''
Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

This file is used to count the reads mapped onto bed format regions
'''

import os,re,argparse,sys,io
import bisect
import gzip,struct
from struct import unpack

plus = re.compile('\+')
minus = re.compile('\-')


def add_region(chrom,outer,inner,regions):
    '''
    Add (unique) start-end in regions
    '''
    if chrom not in regions:
        regions[chrom]={}
    if inner not in regions[chrom]:
        regions[chrom][inner] = set()# [],if do not require unique read
    regions[chrom][inner].add(outer)#append 
    return regions  
    

def get_bed_regions(bedfile,chroms):
    ''' 
    Get tag regions from BED files 
    '''  
     
    infile = open(bedfile,'r')
    # 
    try:
        line = infile.readline()
    except:
        sys.stderr.write('Not a valid BED format of file: {} ! \n\n'.format(bedfile))
        sys.exit(1)
    regions1 = {} # plus strand
    regions2 = {} # minus strand
    while line:
        line = line.strip().split()
        if line[0] in chroms and len(line)>=6:
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            strand = line[5]
            #tag = BED(line[0],line[1],line[2],line[3],line[4],line[5])
            if plus.match(strand) or dot.match(strand):
            #all start(outer) positions with same end(inner) are in a same set/list in plus strand
                regions1 = add_region(chrom,start,end,regions1)
            elif minus.match(strand):
                #all end(outer) positions with same start(inner) are in a same set/list in minus strand
                regions2 = add_region(chrom,end,start,regions2)
        else: #if the line dose not match the bed format
            pass
        line = infile.readline()
    infile.close()
    
    # check if tag info is read into regions  
    if len(regions1)+len(regions2) == 0:
        sys.stderr.write('File <{}> is not of a valid BED format! \n'.format(bedfile))
        sys.exit(1)
    return regions1, regions2



def bam_binary_parse(data):
    '''
    Refer to : https://github.com/taoliu/MACS/blob/master/MACS2/IO/Parser.pyx 
    
    The bitwise flag is made like this:
    dec meaning
    --- -------
    1   paired read
    2   proper pair
    4   query unmapped
    8   mate unmapped
    16  strand of the query (1 -> reverse)
    32  strand of the mate
    64  first read in pair
    128 second read in pair
    256 alignment is not primary
    512 does not pass quality check
    1024    PCR or optical duplicate
    2048    supplementary alignment
    ''' 
    
    if not data:
        return (-1,-1,-1, -1)       
    thisref = unpack('<i',data[0:4])[0]
    thisstart = unpack('<i',data[4:8])[0]
    thistagsize  = unpack('<i',data[16:20])[0]
    (n_cigar_op,  bwflag ) = unpack( '<HH' , data[ 12:16 ] )
    #print(thisref,thisstart,thistagsize)
    #exit(0)
    if bwflag & 4 or bwflag & 512 or bwflag & 256 or bwflag & 2048:
        return ( -1, -1, -1, -1 )       #unmapped sequence or bad sequence or  secondary or supplementary alignment 
    if bwflag & 1:
        # paired read. We should only keep sequence if the mate is mapped
        # and if this is the left mate, all is within  the flag! 
        if not bwflag & 2:
            return ( -1, -1, -1, -1 )   # not a proper pair
        if bwflag & 8:
            return ( -1, -1, -1, -1 )   # the mate is unmapped
        # From Benjamin Schiller https://github.com/benjschiller
        if bwflag & 128:
            # this is not the first read in a pair
            return ( -1, -1, -1, -1 )
        # end of the patch
    # In case of paired-end we have now skipped all possible "bad" pairs
    # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
    # we can treat it as a single read, so just check the strand and calculate its
    # start position... hope I'm right!
    if bwflag & 16:
        # read mapped to minus strand
        l_read_name = unpack( '<B', data[ 8:9 ] )[ 0 ]
        # need to decipher CIGAR string
        for cigar_code in unpack( '<%dI' % (n_cigar_op) , data[ 32 + l_read_name : 32 + l_read_name + n_cigar_op*4 ] ):
            #print(cigar_code>>4,'****')
            if cigar_code & 15 in [ 0, 2, 3, 7, 8 ]:   # they are CIGAR op M/D/N/=/X
                thisstart += cigar_code >> 4
                thistagsize = cigar_code >> 4
        thisstrand = 1
    else:
        thisstrand = 0

    return (thisref,thisstart,thistagsize,thisstrand)


def get_bam_regions(bamfile,chroms):
    ''' 
    Get tag regions from BAM files
    Refer to : https://github.com/taoliu/MACS/blob/master/MACS2/IO/Parser.pyx    
    File is gzip-compatible and binary.
    '''
    gzipped = True
    #try gzip first 
    gfile = gzip.open(bamfile)
    try:
        gfile.read(10)
    except IOError:
        gzipped = False
    gfile.close()   
    if gzipped:
        # open with gzip.open, then wrap it with BufferedReader
        infile = io.BufferedReader(gzip.open(bamfile,mode='rb'))
    else:
        infile = io.open(bamfile,mode='rb')

    fseek = infile.seek
    fread = infile.read
    ftell = infile.tell
    # check the first 3 bytes of BAM file    
    if fread(3).decode('utf-8') == 'BAM':  
        fseek(0)
    else:      
        sys.stderr.write('Not a valid BAM format of file: {} ! \n\n'.format(bamfile))
        sys.exit(1)
    # move to pos 4    
    fseek(4)
    header_len = unpack('<i',fread(4))[0]    
    fseek(header_len + ftell())  
    references = [] 
    # get the number of chromosome
    nc = unpack('<i',fread(4))[0]
    for x in range(nc):   
        #read each chromosome
        nlength = unpack('<i',fread(4))[0]
        refname = fread(nlength)[:-1].decode('utf-8')
        #print(x,refname)
        references.append(refname)
        
        # a: jump over chromosome size?
        fseek(ftell()+4)
        
        # b: or read the chromosome size ? of same fun of forward 4 pos? 
        #len_refname = unpack('<i',fread(4))[0]  
        #print(ftell(),len_refname) 
    #print(references)
    regions1 = {} # plus strand
    regions2 = {} # minus strand    
    while True:
        try:
            entrylength = unpack('<i',fread(4))[0]    
        except struct.error:
            break        

        (chromref,tagpos,tagsize,strand) = bam_binary_parse(fread(entrylength))
        #print(chromref,tagpos,tagsize,strand)
        
        if references[chromref] in chroms:
        #strand -- 0: plus, 1: minus
            if tagpos >= 0 and strand == 0:           
            #all start(outer) positions with same end(inner) are in a same set/list in plus strand
                regions1 = add_region(references[chromref],tagpos,tagpos+tagsize,regions1)
            elif tagpos >= 0 and strand == 1:
            #all end(outer) positions with same start(inner) are in a same set/list in minus strand
                regions2 = add_region(references[chromref],tagpos,tagpos-tagsize,regions2)                
    infile.close()
    
    # check if tag info is read into regions  
    if len(regions1)+len(regions2) == 0:
        sys.stderr.write('File <{}> is not of a valid BAM format! \n'.format(bamfile)) 
        sys.exit(1)  

    return regions1, regions2   


def get_tag_regions(species,format,infile):
    '''
    Get all the unique start-end regions in input bed file,
    separated by chrom and strand
    '''
    chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY', 'chrM']
    
    if format.lower()=="bed":
        return get_bed_regions(infile,chroms)
    elif format.lower()=="bam":
        return get_bam_regions(infile,chroms)


def is_list_sorted(mylist):
    '''
    Check if list is sorted
    '''        
    for i in range(len(mylist)-1):
        if mylist[i] > mylist[i+1]:
            return 0
    return 1   
        	

def get_read_positions(positions,regions,val,fragment_size):
    '''
    Return all the shifted positions of reads 
    '''
    # default fragment_size = 150
    if fragment_size >= 0:
        shift = int(round(fragment_size/2))	
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                positions[chrom].extend([outer+shift*val for outer in regions[chrom][inner]])
    else:
        for chrom in regions.keys():
            if chrom not in positions:
                positions[chrom]=[]
            for inner in regions[chrom].keys():
                for outer in regions[chrom][inner]:
                    shift = int(abs(outer-inner)/2)	
                    positions[chrom].append(outer+shift*val)
    return positions   
	        	    
		
def get_count_on_mapID(start,end,positions,mid,expand):
    '''
    Count the tags/positions on DHS
    '''
    #if is_list_sorted(positions)==0:
    #    positions.sort()
    if start < end and mid:
        middle = int((start+end)*0.5)
        s = bisect.bisect_left(positions,middle-expand)
        e = bisect.bisect_right(positions,middle+expand)
        return e-s, expand*2
        
    elif start < end:
        s = bisect.bisect_left(positions,start-expand)
        e = bisect.bisect_right(positions,end+expand)
        return e-s, end-start+expand*2
            
    else:
        return 0,1


def get_tag_regions(species,format,infile):
    '''
    Get all the unique start-end regions in input bed file,
    separated by chrom and strand
    '''
    chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
         'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
         'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
    
    if format.lower()=="bed":
        return get_bed_regions(infile,chroms)
    elif format.lower()=="bam":
        return get_bam_regions(infile,chroms)


def read_count_on_mapfile(infile,readsfile,species,format,fragmentsize=147,mid=False,expand=0):
    '''
    Count the num of (unique) reads in user-input bed/bam file
    on each of the UDHS -- bart profile
    '''
    #specify the species
    # get the start-end regions of each read in each chrom (as key) and separate by strand

    regions1, regions2 = get_tag_regions(species,format,readsfile)

    # get the counting position of each read(tag), chrom as key
    positions = get_read_positions({},regions1,1,fragmentsize)
    positions = get_read_positions(positions,regions2,-1,fragmentsize)
    total=0
    for chrom in positions:
        positions[chrom].sort() 
        total+=len(positions[chrom])
    # check if positions is NULL
    #if total ==0:
        #sys.stderr.write('Can not read the input bed/bam file!\n')
    # count reads on each DHS

    
    mapfile = open(infile,'r')
    line = mapfile.readline()
    counting = {}
    while line:
        line = line.strip().split()
        #dhs = BED(line[0],line[1],line[2],line[3],line[4],line[5])
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        map_id = line[3]
        if chrom in positions:
            #assert DHS_id not in counting
            nums,region_width = get_count_on_mapID(start,end,positions[chrom],mid,expand)
            counting[map_id]= round(nums*1000000000/(total*region_width),3)
        line = mapfile.readline()
    mapfile.close()
    return counting

def main(infile,readsfile,outfile,species,format,fragmentsize,mid,expand):

    counting = read_count_on_mapfile(infile,readsfile,species,format,fragmentsize,mid,expand)
    with open(outfile,'w') as outf:
        for i in counting:
            outf.write('{}\t{}\n'.format(i,counting[i]))
    
    
    
    
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of bed format, with start position sorted', metavar = '<file>')
    parser.add_argument('-b', '--readsfile', action = 'store', type = str,dest = 'readsfile', help = 'ChIP-seq reads in bed/bam format', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile to write the the RPKM for each region', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
   
    parser.add_argument('-f', '--format', action = 'store', type = str,dest = 'format', help = 'format of file, BED or BAM', metavar = '<str>',required=True)
    parser.add_argument('-g', '--fragmentsize', action = 'store', type = int,dest = 'fragmentsize', help = 'fragmentsize for the shift of reads. Default: 147', metavar = '<int>',default=147)
    parser.add_argument('-m', '--mid', action = 'store_true', dest = 'mid', help = 'whether to use middle side for expansion. Default: False',default=False)
    parser.add_argument('-e', '--expand', action = 'store', type = int,dest = 'expand', help = 'expand of regions. Default: 0', metavar = '<int>',default=0)
    

    args = parser.parse_args()
    if(len(sys.argv))<9:
        parser.print_help()
        sys.exit(1)
  
    main(args.infile,args.readsfile,args.outfile,args.species,args.format,args.fragmentsize,args.mid,args.expand)
