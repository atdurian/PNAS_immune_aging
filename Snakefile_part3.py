import os
import glob

##### Paths
R_PATH=                 '/local10G/resources/R-3.1.2'
SCRIPTS=                '/datastore/debourcy/snakemake_scripts/igh_repertoires'
DATA=                   '/datastore/debourcy/Antibodies/Ellison'
WDIR=                   DATA+'/Snakemake_wdir'
LOGS=                   WDIR+'/logs'
      
##### Sample lists
SAMPLE_LIST = [os.path.basename(FOLDER) for FOLDER in glob.glob(DATA+'/01*')]
PATIENTS = list(set([SAMPLE.split("_")[0] for SAMPLE in SAMPLE_LIST]))
SAMPLES = {PATIENT:[os.path.basename(s) for s in glob.glob(DATA+'/'+PATIENT+'*')] for PATIENT in PATIENTS}
def patient_from_sample(sample):
    for patient in SAMPLES.keys():
      if sample in SAMPLES[patient]:
        return patient
        
RUBL_Nseq = [str(int(1e3))]
GERMLINE_FIELD = ['GERMLINE_IMGT_D_MASK','GERMLINE_IMGT']


##### Rules
rule all:
    input:
        expand(WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab', patient=PATIENTS),
        expand(WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab', sample=SAMPLE_LIST)#,
        expand(WDIR+'/RUBL.Nseq-{N}.{germ_field}.tab', N=RUBL_Nseq, germ_field=GERMLINE_FIELD)

rule define_lineages:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass.tab'
    params: name='define_lineages', partition='long', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{patient}_define_lineages.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/define_lineages.R {input_scratch} &&\
                cp {output_scratch} {WDIR}/")
                
rule compute_RUBL:
    input:
        expand(WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab', patient=PATIENTS),
        expand(WDIR+'/{patient}_personal_hIGHV.fasta', patient=PATIENTS)
    output:
        WDIR+'/RUBL.Nseq-{N}.{germ_field}.tab'
    params: name='compute_RUBL', partition='unrestricted', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/compute_RUBL_Nseq-{N}_{germ_field}.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/compute_RUBL.R {wildcards.N} {wildcards.germ_field} {input_scratch_concat} &&\
               cp {output_scratch} {WDIR}/")     
               
rule count_mutations:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass.tab'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected.tab'
    params: name='count_mutations', partition='long', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{patient}_count_mutations.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/count_mutations_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     

rule count_mutations_nonfunctional:
    input:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass.tab'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected.tab'
    params: name='count_mutations_nonfunctional', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{sample}_count_mutations_nonfunctional.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/count_mutations_nonfunctional_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule remove_error_clouds:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected.tab'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab'
    params: name='remove_error_clouds', partition='long', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{patient}_remove_error_clouds.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/remove_error_clouds.R {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     

rule remove_error_clouds_nonfunctional:
    input:
        lambda wildcards: [WDIR+'/'+wildcards.sample+'_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected.tab']+
                          [WDIR+'/'+patient_from_sample(wildcards.sample)+'_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected.tab']       
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected_remove-error-clouds-pass.tab'
    params: name='remove_error_clouds_nonfunctional', partition='long', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{sample}_remove_error_clouds_nonfunctional.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/remove_error_clouds_nonfunctional.R {input_scratch[0]} {input_scratch[1]} &&\
                cp {output_scratch} {WDIR}/")     
                
rule attach_AA_sequences:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass.tab'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab'
    params: name='attach_AA_sequences', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{patient}_attach_AA_sequences.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/attach_AA_sequences_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule attach_AA_sequences_nonfunctional:
    input:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected_remove-error-clouds-pass.tab'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab'
    params: name='attach_AA_sequences_nonfunctional', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{sample}_attach_AA_sequences_nonfunctional.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/attach_AA_sequences_nonfunctional_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
rule attach_AA_mutation_types:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected.tab'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass_lineage-pass_mutations-pass-corrected_remove-error-clouds-pass_attach-AA-sequences-pass-corrected_attach-AA-mutationtypes-pass-corr.tab'
    params: name='attach_AA_mutation_types', partition='general', cpus_per_task='6', mem_per_cpu='5300', logfile=LOGS+'/{patient}_attach_AA_mutation_types.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/attach_AA_mutation_types_corrected.R {WDIR} {input_scratch} &&\
                cp {output_scratch} {WDIR}/")     
                
                
                