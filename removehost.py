#!/usr/bin/env python

import sys
import os
import glob
import argparse 

def main():
    cli = argparse.ArgumentParser()
    
    cli.add_argument('-i', '--InputFolder', help="Folder containing basecalled Nanopore fastq files. Only files ending in .fastq will be used", required=True)
    cli.add_argument('-r', '--Reference', help="Host Reference fasta or fasta.gz file", required=True)
    cli.add_argument('-o', '--OutputFolder', help="Output Folder. Default is ~/dehost_output/test", required=False, default='~/dehost_output/test')

    cli.add_argument('-t', '--threads', help="Number of threads. Default is 4. More is faster if your computer supports it", type=int, required=False, default=4)
    method = cli.add_mutually_exclusive_group()
    method.add_argument('--Nanopore', action='store_const', dest='seq_method', const='map-ont', default='map-ont')
    method.add_argument('--PacBio', action='store_const', dest='seq_method', const='map-pb')
    args = cli.parse_args()

    files = glob.glob(args.InputFolder+"/*.fastq")

    for i in files:
        base = os.path.splitext(os.path.basename(i))[0]
        #print(base)
        os.system(f"mkdir -p {args.OutputFolder}")
        minimap2_cmd = f"minimap2 -ax {args.seq_method} {args.Reference} {i} -t {args.threads} > {args.OutputFolder}/{base}.sam"
        print(minimap2_cmd)
        os.system(minimap2_cmd)
        samtools_cmd1 = f"samtools view -u -f 4 {args.OutputFolder}/{base}.sam > {args.OutputFolder}/{base}_filtered.sam"
        print(samtools_cmd1)
        os.system(samtools_cmd1)
        samtools_cmd2 = f"samtools bam2fq {args.OutputFolder}/{base}_filtered.sam > {args.OutputFolder}/{base}_filtered.fastq"
        print(samtools_cmd2)
        os.system(samtools_cmd2)
        delete_cmd1 = f"rm {args.OutputFolder}/{base}.sam"
        os.system(delete_cmd1)
        delete_cmd2 = f"rm {args.OutputFolder}/{base}_filtered.sam"
        os.system(delete_cmd2)

if __name__ == "__main__":
    sys.exit(main())
