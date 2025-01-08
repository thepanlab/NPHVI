# A pipline to analyze nanopore direct RNA sequencing data
The human genome and transcriptome reference used Gencode_v45  
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

Calculate the depth/coverage  
```
samtools view -b -F 4 Output_minimap.sam > aligned_minimap.bam
```
```
samtools sort -o aligned_minimap_sorted.bam aligned_minimap.bam
```
```
samtools index aligned_minimap_sorted.bam
```
```
samtools depth aligned_minimap_sorted.bam > coverage.txt
```
Calculate the sum of coverage  
```
awk '{split($1, arr, "|"); ref_name = arr[1]; coverage[ref_name] += $3} END {for (ref in coverage) print ref, coverage[ref]}' coverage.txt > Sum_coverage.txt
```
Calculate the average coverage  
```
awk '{sum+=$3} END {print "Average Coverage = ", sum/NR}' coverage.txt > Avg_coverage.txt
```
Calculate the Breadth of Coverage
```
samtools depth aligned_minimap_sorted.bam | awk '{if($3>0) covered++} END {print covered/1000*100 "%"}' > Breadth_coverage.txt
```
## Calculate transcripts abundance (TPM)
Software: samtools
```
samtools idxstats aligned_minimap_sorted.bam > Count.txt
```
```
awk 'BEGIN {FS="\t"; OFS="\t"} {if($3>0) print $0, ($3/$2)*1000000}' Count.txt > TPM.txt
```
## Segmentation
Software: Nanopolish
https://github.com/jts/nanopolish  
```
module load Eigen/3.2.9-intel-2016a
module load HTSlib/1.9-foss-2018b
export HDF5_PLUGIN_PATH=/ourdisk/hpc/rnafold/dywang/dont_archive/Software/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
```
Index
```
nanopolish index -d directory_fast5 -s sequencing_summary.txt fastq_file.fastq
```
Segmentation  
```
nanopolish eventalign --scale-events --signal-index --threads 50
--reads fastq_file.fastq
--bam aligned_minimap_sorted.bam
--genome gencode.v45.transcripts.fa
--summary sequencing_summary.txt
> event.txt
```
## polyA length calculation
Software: Nanopolish  
```
nanopolish polya -t8 -r fastq_file.fastq -b aligned_minimap_sorted.bam -g gencode.v45.transcripts.fa > polyA_length.txt
```
Selected the results with"pass"
```
grep 'PASS' polyA_length.txt > passed_polyA.txt
```
## m6A modification identification
Software: m6Anet  
https://github.com/GoekeLab/m6anet
```
module load Python/3.8.6-GCCcore-10.2.0
m6anet-dataprep --eventalign event.txt --out_dir m6Anet --n_processes 30
m6anet-run_inference --input_dir m6Anet --out_dir m6Anet --infer_mod_rate --n_processes 0
```
The script Possibility_selection.py was used for selecting the modification that the possibility is higher than 95%  
## Compare the m6A modification  
Software: Xpore  
https://github.com/GoekeLab/xpore  
```
xpore dataprep --eventalign event.txt --out_dir Xpore
```
```
xpore diffmod --config config.yml
```
























