# A pipline to analyze nanopore direct RNA sequencing data
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
