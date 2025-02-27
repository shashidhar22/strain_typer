import os
import sys
import csv
import glob
import shutil
import logging
import argparse
import tempfile
import subprocess
import pandas
from Bio import SeqIO
global relpath
relpath = os.path.abspath(os.path.dirname(__file__))

def gff_parser(gfffile,genelist,outdir):
    with open(genelist) as gen_list:
        genes = [vals.strip() for vals in gen_list.readlines()]
    sero_group  = {'csb':'Nm_serob','csc':'Nm_seroC','csy':'Nm_seroY','csw':'Nm_seroW','csxB':'Nm_seroX','csaB':'Nm_seroA','siaD':'Nm_serob'}
    outfile = csv.writer(open(outdir+'/moltyping.csv','w'),delimiter='\t')
    outfile.writerow(['Gene','Status','Serogroup/Serotype'])
    for gene_name in genes:
        if gene_name in open(gfffile).read():
              if gene_name in sero_group.keys():
                  outfile.writerow([gene_name,'Present',sero_group[gene_name]])
              else:
                  outfile.writerow([gene_name,'Present','NA'])
        else:
              outfile.writerow([gene_name,'Absent','NA'])
    return

def mol_typer(gene_file, inpfile, outdir):
    """Module to perform molecular characterization based given set of genes"""
    #create gene index
    fasta_file = SeqIO.parse(open(gene_file),'fasta')
    gene_index = dict()
    mask = (1<<(11*2))-1
    exist = str()
    for count,sequences in enumerate(fasta_file):
        header = sequences.id
        sequence = str(sequences.seq)
        kmer = 0
        bit_counter = 0
        for nuc in sequence:
            kmer =  kmer << 2         #left shift k-kmer 
            if nuc == 'A':            #add the value of the character using bitwise OR
                kmer = kmer | 0x00
            elif nuc == 'C':
                kmer = kmer | 0x01
            elif nuc == 'G':
                kmer = kmer | 0x02
            elif nuc == 'T':
                kmer = kmer | 0x03
            else:
                bit_counter = 0
            if bit_counter == 11:   #if length equals k-mer length, store k-mer
                try:
                    exist = gene_index[kmer & mask]  #k-mer to transcript index mapping
                    del gene_index[kmer & mask]
                except KeyError:
                    gene_index[int(kmer & mask)] = header
            else:
                bit_counter += 1
    file_type = detect_type(inpfile)
    run_kan = subprocess.Popen(['count','-d','1','-l','1','-m','dec','--countfilter=kmercount:c>1','-k','11','-f',file_type,'-o',outdir+'/tmp.kc',inpfile], stdout=subprocess.PIPE, shell=False)
    run_kan.wait()
    gene_evidence = dict()
    for lines in csv.reader(open(outdir+'/tmp.kc'),delimiter='\t'):
        if int(lines[0]) in gene_index.keys():
            try:
                gene_evidence[gene_index[int(lines[0])]] += 1
            except KeyError:
                gene_evidence[gene_index[int(lines[0])]]  = 1
    print(len(set(gene_index.values())))
    print(len(gene_evidence))
    os.remove(outdir+'/tmp.kc')
    return

        
def run_mlst(species, inpfile, outdir):
    """Module to run MLST analysis on input sequences for NM and HI species"""
    MLST = {'NM':[relpath+'/MLST/Neiserria_meningitidis/Profile/neisseria.txt',relpath+'/MLST/Neiserria_meningitidis/Sequences'],
            'HI':[relpath+'/MLST/Haemophillus_influenzae/Profile/hinfluenzae.txt',relpath+'/MLST/Haemophillus_influenzae/Sequences']}
    temp_folder = outdir+'/tmp_fasta'
    os.mkdir(temp_folder)
    for files in inpfile:
        shutil.copy(files,temp_folder)
    run_mlst = subprocess.Popen([relpath+'/run_MLST.py','-i',MLST[species][1],'-g',temp_folder,'-p',MLST[species][0],'-o',outdir+'/MLST'], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    run_mlst.wait()
    shutil.copy(outdir+'/MLST/MLST.tab',outdir+'/'+species+'_MLST.tab')
    for files in glob.glob(temp_folder+'/*'):
        os.remove(files)
    os.rmdir(temp_folder)
    for files in glob.glob(outdir+'/MLST/*'):
        os.remove(files)
    os.rmdir(outdir+'/MLST')
    run_mlst.wait()
    return

def detect_type(inpfile):
    """Detect input file type"""
    seqfile = open(inpfile)
    seqline = seqfile.read()
    seqfile.close()
    if ">" == seqline[0]:
        return('fasta')
    else:
        return('fastq')
 
def compare(indir,outdir,loglevel,ref,mode):
    """Species delineation and MLST analysis"""
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
         raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level, format='%(levelname)s:%(asctime)s:%(message)s', datefmt='%m/%d/%Y;%I:%M:%S')
    logging.info('Starting strain typing')
    #Gathering file paths and determining output file names and path
    if ref == None:
        logging.error('Please provide reference files')
        logging.error('Exiting program')
        sys.exit()
    refiles = glob.glob(os.path.abspath(ref)+'/*kc')
    species = {os.path.basename(fasta): os.path.splitext(os.path.basename(fasta))[0].split('_')[0] for fasta in refiles}
    inpbase = [os.path.basename(fasta) for fasta in indir]
    basenames = [os.path.splitext(os.path.basename(fasta))[0] for fasta in indir]
    kcbase = [vals+'.kc' for vals in basenames]
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    klen = '11'
    for files in refiles:
        shutil.copy(files,outdir)
    logging.info('k-merizing '+str(len(basenames))+' files')
    #Running KAnalyze on all input files; can be parallelized
    filterval = {'fasta':'1','fastq':'3'}
    for fasta,base in zip(indir,basenames):
        file_type = detect_type(fasta)
        run_kan = subprocess.Popen([relpath+'/kanalyze-0.9.7/count','-d','1','-l','1','--countfilter=kmercount:c>'+filterval[file_type],'-k',klen,'-f',file_type,'-o',outdir+'/'+base+'.kc',fasta], stdout=subprocess.PIPE, shell=False)
        run_kan.wait()
    logging.info('Running Rspecies')
    #Running Rspecies
    run_rspecies = subprocess.Popen(['Rscript',relpath+'/Rspecies.R',klen,outdir,outdir+'/'],shell=False)
    run_rspecies.wait()
    logging.info('Tree created')
    phylo_tree = pandas.read_table(outdir+'/dist.csv',header=0,index_col=0,sep=',',skipinitialspace=True)
    species_file = csv.writer(open(outdir+'/species.csv','w'),delimiter='\t')
    species_file.writerow(['Assembly','Closest Reference','Predicted Species','Distance to Reference'])
    mlst_classifier = {'NM':[],'HI':[]}
    logging.info('Determining species')
    #Identifying minimum distance reference
    for values,inpfile in zip(kcbase,indir):
        minimum =1
        minval = int()
        for vals in species.keys():
            if phylo_tree[values][vals] < minimum: 
                minval = vals
                minimum = phylo_tree[values][vals]
        species_file.writerow([values,minval,species[minval],minimum ])
        try:
            mlst_classifier[species[minval]].append(inpfile)
        except KeyError:
            continue           
    logging.info('Clearing temporary files')
    #Deleting temporary files
    for kcfiles in glob.glob(outdir+'/*.kc'):
        os.remove(kcfiles)
    #Running MLST analsyis
    if mode == 'both' and file_type == 'fasta':
        logging.info('Running MLST analysis')
        for species in mlst_classifier:
            run_mlst(species,mlst_classifier[species], outdir)
    else:
        logging.error('Can\'t run MLST on fastq files')
        logging.error('Exiting program.')
        sys.exit()
    logging.info('Analysis done')
    return
      

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='strain_typer')
    parser.add_argument('-i','--input_file',dest='indir',nargs="+",type=str, help='Path to files to be analyzed')
    parser.add_argument('-r','--ref_dir',dest='refdir',type=str,help='Path to refrence kc files')
    parser.add_argument('-o','--output_dir',dest='outdir',type=str, help='Output directory path')
    parser.add_argument('-g','--genefile',dest='gene',type=str,help='Molecular typing gene file')
    parser.add_argument('-m','--mode',dest='mode',type=str,help='Mode of analysis',choices=['species','mlst','both','anno'],default='both')
    parser.add_argument('-s','--species',dest='species',type=str,help='Species against which MLST needs to be performed. Must be specified in MLST mode')
    parser.add_argument('-l','--log',type=str,default='info',choices=['info','debug','error'],help='Verbosity parameter')
    parser.add_argument('-a','--annotation',dest='gff',type=str,help='GFF file')
    parser.add_argument('-e','--genelist',dest='gl',type=str,help='Gene list file')
    args = parser.parse_args()
    if args.mode == 'both' or args.mode == 'species':
        compare(args.indir, args.outdir,args.log,args.refdir,args.mode)
    elif args.mode == 'mlst' and args.species != None:
        run_mlst(args.species,args.indir,args.outdir)
    elif args.mode == 'mlst' and args.species == None:
        print('Species not specified for MLST analysis. Terminating strain typer.')
        sys.exit()
    elif args.mode == 'anno':
        gff_parser(args.gff,args.gl,args.outdir)
#    mol_typer(args.gene,args.indir[0],args.outdir)
