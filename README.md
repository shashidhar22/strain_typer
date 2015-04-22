# Strain typing tool
#### version-0.9.1
This tool uses an alignment free methodlogy to identify the species of the input fasta or fastq file.
Currently species identification is limited to *Neisseria meningitidis*,*Haemophilus influenzae* and 
*Haemophilus haemolyticus*.
The tool also identifies MLST using the [run_MLST](https://github.com/widdowquinn/scripts) script from
Leighton Pritchard's bioinformatics toolkit.

##Prerequisites:
1. R
2. Python (preferably python3.4)
3. Blast
4. Java

##Software and Packages:
1. [KAnalyze](http://sourceforge.net/projects/kanalyze/)
2. [run_MLST](https://github.com/widdowquinn/scripts)
3. kmerDistance R package

##Execution:
To run strain typer, the strain_typer.py needs to be executed.
```python
python3.4 strain_typer.py -h
usage: strain_typer [-h] [-i INDIR [INDIR ...]] [-r REFDIR] [-o OUTDIR]
                    [-g GENE] [-m {species,mlst,both,anno}] [-s SPECIES]
                    [-l {info,debug,error}] [-a GFF] [-e GL]

optional arguments:
  -h, --help            show this help message and exit
  -i INDIR [INDIR ...], --input_file INDIR [INDIR ...]
                        Path to files to be analyzed
  -r REFDIR, --ref_dir REFDIR
                        Path to refrence kc files
  -o OUTDIR, --output_dir OUTDIR
                        Output directory path
  -g GENE, --genefile GENE
                        Molecular typing gene file
  -m {species,mlst,both}, --mode {species,mlst,both}
                        Mode of analysis
  -s SPECIES, --species SPECIES
                        Species against which MLST needs to be performed. Must
                        be specified in MLST mode
  -l {info,debug,error}, --log {info,debug,error}
                        Verbosity parameter
  -a GFF, --annotation GFF
                        GFF file
  -e GL, --genelist GL  Gene list file
```

##Examples:
```python
time python3.4 strain_typer.py -i ./SPADES/corrected_fasta/M08706.fa -o ./pipeline_test -r ref_kc/
```
