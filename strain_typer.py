import os
import sys
import csv
import glob
import logging
import argparse
import subprocess
import pandas

def main(indir,ftype, outdir,loglevel,ref):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
         raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level, format='%(levelname)s:%(asctime)s:%(message)s', datefmt='%m/%d/%Y;%I:%M:%S')
    logging.info('Starting strain typing')
    files = glob.glob(os.path.abspath(ref+'/*fna'))
    refiles = glob.glob(os.path.abspath(ref+'/*fna'))
    species = {os.path.splitext(os.path.basename(fasta))[0]: os.path.splitext(os.path.basename(fasta))[0].split('_')[0] for fasta in refiles}
    files += indir
    inpbase = [os.path.splitext(os.path.basename(fasta))[0] for fasta in indir]
    basenames = [os.path.splitext(os.path.basename(fasta))[0] for fasta in files]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    klen = '11'
    logging.info('k-merizing '+str(len(basenames))+' files')
    for fasta,base in zip(files,basenames):
#        print(' '.join(['/home/sravishankar9/ComparativeGenomics/kanalyze-0.9.7/count','-d','6','-l','6','--countfilter=kmercount:c>1','-k',klen,'-f',ftype,'-o',outdir+'/'+base+'.kc',fasta]))
        if fasta in refiles or ftype == 'fasta':
            run_kan = subprocess.Popen(['count','-d','1','-l','1','--countfilter=kmercount:c>1','-k',klen,'-f','fasta','-o',outdir+'/'+base+'.kc',fasta], stdout=subprocess.PIPE, shell=False)
            run_kan.wait()
        elif ftype == 'fastq':
            run_kan = subprocess.Popen(['count','-d','1','-l','1','--countfilter=kmercount:c>1','-k',klen,'-f','fastq','-o',outdir+'/'+base+'.kc',fasta], stdout=subprocess.PIPE, shell=False)
            run_kan.wait()
        else:
            logging.error('Illegal file type')
            logging.error('Exiting program')
            os.exit()
    logging.info('Running Rspecies')
    run_rspecies = subprocess.Popen(['Rscript','/home/shashidhar/code_repos/strain_typer/Rspecies.R',klen,outdir,outdir+'/'],shell=False)
    run_rspecies.wait()
    logging.info('Tree created')
    phylo_tree = pandas.read_table(outdir+'/dist.csv',header=0,index_col=0,sep=',',skipinitialspace=True)
    species_file = csv.writer(open(outdir+'/species.csv','w'),delimiter='\t')
    species_file.writerow(['Assembly','Closest Reference','Predicted Species','Distance to Reference'])
    for values in inpbase:
        minimum =1
        minval = int()
        for vals in range(len(phylo_tree[values+'.kc'])):
            if phylo_tree[values+'.kc'][vals] < minimum and phylo_tree[values+'.kc'][vals] != 0.0 and 'SPADES' not in phylo_tree.iloc[vals].name:
                minval = vals
                minimum = phylo_tree[values+'.kc'][vals]
        print(values,phylo_tree.iloc[minval].name,species[phylo_tree.iloc[minval].name.split('.')[0]], phylo_tree[values+'.kc'][minval])
        species_file.writerow([values,phylo_tree.iloc[minval].name,species[phylo_tree.iloc[minval].name.split('.')[0]], phylo_tree[values+'.kc'][minval]])
    logging.info('Analysis done')
    return
      

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='strain_typer')
    parser.add_argument('-i','--input_file',dest='indir',nargs="+",type=str, help='Path to files to be analyzed')
    parser.add_argument('-r','--ref_dir',dest='refdir',type=str,help='Path to refrence fasta files')
    parser.add_argument('-t','--type',dest='ftype',type=str, help='File type[fasta/fastq]')
    parser.add_argument('-o','--output_dir',dest='outdir',type=str, help='Output directory path')
    parser.add_argument('-l','--log',type=str,default='info',choices=['info','debug','error'],help='Verbosity parameter')
    args = parser.parse_args()
    main(args.indir, args.ftype,args.outdir,args.log,args.refdir)
