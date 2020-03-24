#!/usr/bin/env python

import sys
import os
import glob
from gooey import Gooey, GooeyParser

@Gooey(default_size=(720, 750))
        

def main():
    cli = GooeyParser(description="Remove Host Reads from Nanopore Fastq Files")
    required_args = cli.add_argument_group("Input Output Location", gooey_options={'columns': 1})
    required_args.add_argument('--InputFolder', help="Folder containing basecalled Nanopore fastq files. Only files ending in .fastq will be used", required=True, widget='DirChooser')
    required_args.add_argument('--Reference', help="Host Reference fasta or fasta.gz file", required=True, widget='FileChooser')
    required_args.add_argument('--OutputFolder', help="Output Folder", required=False, default='~/dehost_output/test')

    parser = cli.add_argument_group("Optional Arguments", gooey_options={'columns': 1})
    parser.add_argument('--threads', help="Number of threads. More is faster if your computer supports it", type=int, required=False, default=4)

    args = cli.parse_args()

    files = glob.glob(args.InputFolder+"/*.fastq")

    for i in files:
        base = os.path.splitext(os.path.basename(i))[0]
        #print(base)
        os.system(f"mkdir -p {args.OutputFolder}")
        minimap2_cmd = f"minimap2 -ax map-ont {args.Reference} {i} -t {args.threads} > {args.OutputFolder}/{base}.sam"
        print(minimap2_cmd)
        os.system(minimap2_cmd)
        samtools_cmd1 = f"samtools view -u -f 4 {args.OutputFolder}/{base}.sam > {args.OutputFolder}/{base}_filtered.sam"
        print(samtools_cmd1)
        os.system(samtools_cmd1)
        samtools_cmd2 = f"samtools bam2fq {args.OutputFolder}/{base}_filtered.sam > {args.OutputFolder}/{base}_filtered.fastq"
        print(samtools_cmd2)
        os.system(samtools_cmd2)
        delete_cmd = f"rm {args.OutputFolder}/*.sam"
        os.system(delete_cmd)

if __name__ == "__main__":
    sys.exit(main())
