# A pipline to analyze nanopore direct RNA sequencing data
The humance genome and transcriptome reference used Gencode_v45  
## Basecalling
Software: Guppy  
https://timkahlke.github.io/LongRead_tutorials/BS_G.html  
```
guppy_basecaller --compress_fastq -i *.fast5 --save_path .fastQ --config config.cfg
```
Using the following way to run jobs as array
```
   > #SBATCH --output=array_%A_%a.out

   > #SBATCH --error=array_%A_%a.err

   > #SBATCH --array=239-620

   echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
   
   guppy_basecaller -i / _${SLURM_ARRAY_TASK_ID}.fast5 --save_path /fastQ_1 --config config.cfg
```
Output: file.fastq, sequencing_summary.txt  

Merge all the generated fastQ files together
```
cat *.fastq > all.fastq
```
## Alignment
Software: Minimap2  
https://github.com/lh3/minimap2    
```
minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam
```
## Quality control
Software: Marginalign, samtools  
https://github.com/benedictpaten/marginAlign  
https://github.com/samtools/samtools  

Build the virtual environment  
```
 module load Python/2.7.14-intel-2018a
 cd marginAlign
 make clean
 make
 virtualenv --python=python2.7 --no-site-packages --distribute env && source env/bin/activate && pip install -r requirements.txt
```
Before calculating by marginalign, unique the headers in the fastQ file is required (Script was provided in the Script folder)
```
uniquifyFastq input.fastQ input_unique_header
```
To calculate the median/avg/min/max identity of reads in a sam file
```
marginStats
--localAlignment input.sam read.fastq
reference.fasta
--readIdentity
--alignmentIdentity
--readCoverage --mismatchesPerAlignedBase
--deletionsPerReadBase
--insertionsPerReadBase
--readLength
--printValuePerReadAlignment
```
To calculate the substitution  
```
marginAlign/scripts/substitutions _unique.fastq reference.fa _minimap_pass.sam output_directory
```
Using the script "substitution_plot.R" and "identityPlots.R" to visualize the results  




