# PNAS_immune_aging
1) 
* Install snakemake. Then look at  Snakefile_part1.py:
* Put the R-scripts in the SCRIPTS folder
* Put IMGT_hIGHV.fasta, IMGT_hIGHJ.fasta, IMGT_hIGHD.fasta in the folder IMGT_GERMLINES
* In general path names in the Snakefiles as appropriate for you
* Install the R-packages required at the different stages of the workflow.
* Then execute Snakefile_part1.py, e.g.:
snakemake -s Snakefile_part1.py -j 999 --latency-wait 60 --cluster "sbatch --job-name {params.name} --ntasks=1 --cpus-per-task={params.cpus_per_task} --partition={params.partition} --mem-per-cpu={params.mem_per_cpu} --output={params.logfile}"

2) After Snakefile_part1.py has finished, submit the {patient}_minCONSCOUNT2*.fasta files to IMGT/High-V Quest.

3) Download the IMGT/High-V Quest output ".txz" files into the IMGT_RESULTS folder and execute Snakefile_part2.py, e.g.:
snakemake -s Snakefile_part2.py -j 999 --latency-wait 60 --cluster "sbatch --job-name {params.name} --ntasks=1 --cpus-per-task={params.cpus_per_task} --partition={params.partition} --mem-per-cpu={params.mem_per_cpu} --output={params.logfile}"

4) Execute Snakefile_part3.py and Snakefile_part4.py similarly:
snakemake -s Snakefile_part3.py -j 999 --latency-wait 60 --cluster "sbatch --job-name {params.name} --ntasks=1 --cpus-per-task={params.cpus_per_task} --partition={params.partition} --mem-per-cpu={params.mem_per_cpu} --output={params.logfile}"
snakemake -s Snakefile_part4.py -j 999 --latency-wait 60 --cluster "sbatch --job-name {params.name} --ntasks=1 --cpus-per-task={params.cpus_per_task} --partition={params.partition} --mem-per-cpu={params.mem_per_cpu} --output={params.logfile}"

4) Now you should have all the files necessary for running the figure-making scripts in R:
MakeFig1BC.R, MakeFig1D.R, MakeFig2A.R, MakeFig2BCDE.R, MakeFig3ABC.R, MakeFig3D.R, MakeFig3E.R, MakeFig3F.R, MakeFig4AB.R, MakeFig4C.R, MakeFig4D.R
