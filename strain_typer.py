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

def run_mlst(species, inpfile, outdir):
    MLST = {'NM':['/home/shashidhar/code_repos/strain_typer/MLST/Neiserria_meningitidis/Profile/neisseria.txt','/home/shashidhar/code_repos/strain_typer/MLST/Neiserria_meningitidis/Sequences'],
            'HI':['/home/shashidhar/code_repos/strain_typer/MLST/Haemophillus_influenzae/Profile/hinfluenzae.txt','/home/shashidhar/code_repos/strain_typer/MLST/Haemophillus_influenzae/Sequences']}
    temp_folder = outdir+'/tmp_fasta'
    os.mkdir(temp_folder)
    basename = os.path.splitext(os.path.basename(inpfile))[0]
    shutil.copy(inpfile,temp_folder)
    run_mlst = subprocess.Popen(['/home/shashidhar/code_repos/strain_typer/run_MLST.py','-i',MLST[species][1],'-g',temp_folder,'-p',MLST[species][0],'--blast_exe','/usr/bin/blastn','-o',outdir+'/MLST'], shell=False)
    run_mlst.wait()
    shutil.copy(outdir+'/MLST/MLST.tab',outdir+'/'+basename+'_MLST.tab')
    for files in glob.glob(temp_folder+'/*'):
        os.remove(files)
    os.rmdir(temp_folder)
    for files in glob.glob(outdir+'/MLST/*'):
        os.remove(files)
    os.rmdir(outdir+'/MLST')
    run_mlst.wait()
    return

def detect_type(inpfile):
    seqfile = open(inpfile)
    seqline = seqfile.read()
    seqfile.close()
    if ">" in seqline:
        return('fasta')
    else:
        return('fastq')
 
def typer(refpath, inpfile ,outdir):
    #make temp directory
    os.mkdir(outdir)
    ref_files = glob.glob(refpath+'*.kc')
    ref_base = [os.path.splitext(os.path.basename(vals))[0] for vals in ref_files]
    species = {os.path.basename(vals): os.path.splitext(os.path.basename(vals))[0].split('_')[0] for vals in ref_files}
    temp_dir = os.path.abspath(outdir)+'/tmp_kc'
    os.mkdir(temp_dir)
    for files in ref_files:
        shutil.copy(files,temp_dir)
    inpbase = os.path.splitext(os.path.basename(inpfile))[0] 
    file_type = detect_type(inpfile)
    run_kan = subprocess.Popen(['count','-d','1','-l','1','--countfilter=kmercount:c>1','-k','11','-f',file_type,'-o',temp_dir+'/'+inpbase+'.kc',inpfile], stdout=subprocess.PIPE, shell=False)
    run_kan.wait()
    run_species = subprocess.Popen(['Rscript','/home/shashidhar/code_repos/strain_typer/Rspecies.R','11',temp_dir+'/',outdir+'/'],shell=False)
    run_species.wait()
    phylo_tree = pandas.read_table(outdir+'/dist.csv',header=0,index_col=0,sep=',',skipinitialspace=True)
    minimum = 100
    minval = str()
    for vals in range(len(phylo_tree[inpbase+'.kc'])):
        if phylo_tree[inpbase+'.kc'][vals] < minimum and phylo_tree.iloc[vals].name != inpbase+'.kc':
            print(phylo_tree.iloc[vals].name)
            minval = phylo_tree.iloc[vals].name
            minimum = phylo_tree[inpbase+'.kc'][vals]
    for files in glob.glob(temp_dir+'/*'):
        os.remove(files)
    os.rmdir(temp_dir)
    print(inpbase,minval,species[minval], phylo_tree[inpbase+'.kc'][minval])
    run_mlst(species[minval],inpfile,outdir)
    return

def compare(indir,ftype, outdir,loglevel,ref):
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
         raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level, format='%(levelname)s:%(asctime)s:%(message)s', datefmt='%m/%d/%Y;%I:%M:%S')
    logging.info('Starting strain typing')
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
    for fasta,base in zip(indir,basenames):
        file_type = detect_type(fasta)
        run_kan = subprocess.Popen(['count','-d','1','-l','1','--countfilter=kmercount:c>1','-k',klen,'-f',file_type,'-o',outdir+'/'+base+'.kc',fasta], stdout=subprocess.PIPE, shell=False)
        run_kan.wait()
    logging.info('Running Rspecies')
    run_rspecies = subprocess.Popen(['Rscript','/home/shashidhar/code_repos/strain_typer/Rspecies.R',klen,outdir,outdir+'/'],shell=False)
    run_rspecies.wait()
    logging.info('Tree created')
    phylo_tree = pandas.read_table(outdir+'/dist.csv',header=0,index_col=0,sep=',',skipinitialspace=True)
    species_file = csv.writer(open(outdir+'/species.csv','w'),delimiter='\t')
    species_file.writerow(['Assembly','Closest Reference','Predicted Species','Distance to Reference'])
    for values in kcbase:
        minimum =1
        minval = int()
        for vals in species.keys():
            #print(phylo_tree.iloc[vals].name)
            if phylo_tree[values][vals] < minimum: # and phylo_tree.iloc[vals].name not in inpbase:
                print(vals)
                #phylo_tree[values+'.kc'][vals] != 0.0 and 'SPADES' not in phylo_tree.iloc[vals].name:
                minval = vals
                minimum = phylo_tree[values][vals]
        print(values,minval,species[minval], minimum)
        species_file.writerow([values,minval,species[minval],minimum ])
    for kcfiles in glob.glob(outdir+'/*.kc'):
        os.remove(kcfiles)
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
#    compare(args.indir, args.ftype,args.outdir,args.log,args.refdir)
    typer(args.refdir,args.indir[0],args.outdir)
#    run_mlst('NM',args.indir[0],args.outdir)
