import os
import glob

##### Paths and sample lists: Ellison
QIIME_VIRTUALENV=       '/local10G/debourcy/miniconda3/envs/qiime1'
R_PATH=                 '/local10G/resources/R-3.1.2'
SCRIPTS=                '/datastore/debourcy/snakemake_scripts/igh_repertoires'
DATA=                   '/datastore/debourcy/Antibodies/Ellison'
WDIR=                   DATA+'/Snakemake_wdir'
LOGS=                   WDIR+'/logs'
      
SAMPLE_LIST = [os.path.basename(FOLDER) for FOLDER in glob.glob(DATA+'/01*')]
PATIENTS = list(set([SAMPLE.split("_")[0] for SAMPLE in SAMPLE_LIST]))   
CUP_r = [100,20]

RUBL_Nseq = [str(int(1e4))]
GERMLINE_FIELD = ['GERMLINE_IMGT']
SEED = [1,2,3]


##### Rules
rule all:
    input:
        expand(WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_CUPs-{rval}.tab', rval=CUP_r, patient=PATIENTS),
        expand(WDIR+'/RUBL.seed-{seed}.Nseq-{N}.{germ_field}.tab', seed=SEED, N=RUBL_Nseq, germ_field=GERMLINE_FIELD),
        expand(WDIR+'/RUBL2.unnormalized.seed-{seed}.Nseq-{N}.{germ_field}.tab', seed=[1], N=RUBL_Nseq, germ_field=GERMLINE_FIELD),
        expand(WDIR+'/RUBL2.seed-{seed}.Nseq-{N}.{germ_field}.tab', seed=[1], N=RUBL_Nseq, germ_field=GERMLINE_FIELD),


rule compute_CUPs:
    ### input file was made using script "write_biom.R"
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_OTU-table.biom'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_CUPs-{rval}.tab'
    params: name='compute_CUPs', partition='unrestricted', cpus_per_task='1', mem_per_cpu='5300', logfile=LOGS+'/{patient}_compute_CUPs-{rval}.log'
    run:
        input_scratch = os.path.basename(input[0])  
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               set +o nounset &&\
               source {QIIME_VIRTUALENV}/bin/activate {QIIME_VIRTUALENV} &&\
               set -o nounset &&\
               python {QIIME_VIRTUALENV}/bin/conditional_uncovered_probability.py -i {input_scratch} -o {output_scratch} -r {wildcards.rval} -m lladser_pe &&\
               cp {output_scratch} {WDIR}/")
               
rule compute_RUBL_specifiedSeed:
    input:
        expand(WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab', patient=PATIENTS),
        expand(WDIR+'/{patient}_personal_hIGHV.fasta', patient=PATIENTS)
    output:
        WDIR+'/RUBL.seed-{seed}.Nseq-{N}.{germ_field}.tab'
    params: name='compute_RUBL_specifiedSeed', partition='unrestricted', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/compute_RUBL_specifiedSeed_seed-{seed}_Nseq-{N}_{germ_field}.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/compute_RUBL_specifiedSeed.R {wildcards.seed} {wildcards.N} {wildcards.germ_field} {input_scratch_concat} &&\
               cp {output_scratch} {WDIR}/")     
               
rule compute_RUBL2_unnormalized:
    input:
        expand(WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab', patient=PATIENTS),
        expand(WDIR+'/{patient}_personal_hIGHV.fasta', patient=PATIENTS)
    output:
        WDIR+'/RUBL2.unnormalized.seed-{seed}.Nseq-{N}.{germ_field}.tab',
        WDIR+'/RUBL2.seed-{seed}.Nseq-{N}.{germ_field}.tab'
    params: name='compute_RUBL2_unnormalized', partition='unrestricted', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/compute_RUBL2_unnormalized_seed-{seed}_Nseq-{N}_{germ_field}.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/compute_RUBL2_unnormalized.R {wildcards.seed} {wildcards.N} {wildcards.germ_field} {input_scratch_concat} &&\
               cp {output_scratch_concat} {WDIR}/")       
               
              

               
                
