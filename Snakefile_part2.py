import os
import glob

##### Paths
VIRTUALENV=             '/local10G/debourcy/miniconda3/envs/snakes'
CHANGEO_PATH=           '/local10G/debourcy/tools/changeo-0.3.3/bin'
R_PATH=                 '/local10G/resources/R-3.1.2'
IMGT_GERMLINES=         '/local10G/debourcy/tools/tigger_germlines'
SCRIPTS=                '/datastore/debourcy/snakemake_scripts/igh_repertoires'
DATA=                   '/datastore/debourcy/Antibodies/Ellison'
IMGT_RESULTS=           DATA+'/IMGT_output' 
WDIR=                   DATA+'/Snakemake_wdir'
LOGS=                   WDIR+'/logs'
    
##### Sample lists
SAMPLE_LIST = [os.path.basename(FOLDER) for FOLDER in glob.glob(DATA+'/01*')]
PATIENTS = list(set([SAMPLE.split("_")[0] for SAMPLE in SAMPLE_LIST]))
SAMPLES = {PATIENT:[os.path.basename(s) for s in glob.glob(DATA+'/'+PATIENT+'*')] for PATIENT in PATIENTS}
IMGT_PIECES = {SAMPLE:[os.path.splitext(os.path.basename(s))[0] for s in glob.glob(IMGT_RESULTS+'/'+SAMPLE+'*.txz')] for SAMPLE in SAMPLE_LIST}
def patient_from_sample(sample):
    for patient in SAMPLES.keys():
      if sample in SAMPLES[patient]:
        return patient

##### Rules
rule all:
    input:
        expand(WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab', patient=PATIENTS),
        expand(WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass.tab', sample=SAMPLE_LIST)
    
rule create_germlines_nomask_nonfunctional:
    input:
        lambda wildcards: [WDIR+'/'+wildcards.sample+'_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass.tab']+
                          [WDIR+'/'+patient_from_sample(wildcards.sample)+'_personal_hIGHV.fasta']              
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-pass.tab',
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass_nomask_germ-fail.tab'
    params: name='create_germlines_nomask_nonfunctional', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_create_germlines_nomask_nonfunctional.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        prefix = os.path.splitext(input_scratch[0])[0]+'_nomask'
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir germline_repo &&\
               cp {input_scratch[1]} germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHD.fasta germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHJ.fasta germline_repo &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/CreateGermlines.py -d {input_scratch[0]} --failed -r germline_repo -g full --outname {prefix} --vf V_CALL_GENOTYPED --sf SEQUENCE_IMGT &&\
               cp {output_scratch_concat} {WDIR}/")

rule create_germlines_nomask:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass.tab',
        WDIR+'/{patient}_personal_hIGHV.fasta'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-pass.tab',
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass_nomask_germ-fail.tab'
    params: name='create_germlines_nomask', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{patient}_create_germlines_nomask.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        prefix = os.path.splitext(input_scratch[0])[0]+'_nomask'
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir germline_repo &&\
               cp {input_scratch[1]} germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHD.fasta germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHJ.fasta germline_repo &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/CreateGermlines.py -d {input_scratch[0]} --failed -r germline_repo -g full --outname {prefix} --vf V_CALL_GENOTYPED --sf SEQUENCE_IMGT &&\
               cp {output_scratch_concat} {WDIR}/")

rule create_germlines_nonfunctional:
    input:
        lambda wildcards: [WDIR+'/'+wildcards.sample+'_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass.tab']+
                          [WDIR+'/'+patient_from_sample(wildcards.sample)+'_personal_hIGHV.fasta']      
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-pass.tab',
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass_germ-fail.tab'
    params: name='create_germlines', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_create_germlines_nonfunctional.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir germline_repo &&\
               cp {input_scratch[1]} germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHD.fasta germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHJ.fasta germline_repo &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/CreateGermlines.py -d {input_scratch[0]} --failed -r germline_repo --vf V_CALL_GENOTYPED --sf SEQUENCE_IMGT &&\
               cp {output_scratch_concat} {WDIR}/")
               
rule create_germlines:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass.tab',
        WDIR+'/{patient}_personal_hIGHV.fasta'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-pass.tab',
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass_germ-fail.tab'
    params: name='create_germlines', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{patient}_create_germlines.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir germline_repo &&\
               cp {input_scratch[1]} germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHD.fasta germline_repo &&\
               cp {IMGT_GERMLINES}/IMGT_hIGHJ.fasta germline_repo &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/CreateGermlines.py -d {input_scratch[0]} --failed -r germline_repo --vf V_CALL_GENOTYPED --sf SEQUENCE_IMGT &&\
               cp {output_scratch_concat} {WDIR}/")

rule correct_nonfunctional_Vcalls:
    input:
        lambda wildcards: [WDIR+'/'+wildcards.sample+'_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE.tab']+
                          [WDIR+'/'+patient_from_sample(wildcards.sample)+'_personal_hIGHV.fasta']
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE_V-CALL-GENOTYPED-pass.tab'
    params: name='correct_nonfunctional_Vcalls', partition='general', nproc='4', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_correct_nonfunctional_Vcalls.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/correct_nonfunctional_Vcalls.R {input_scratch[0]} {input_scratch[1]} &&\
               cp {output_scratch_concat} {WDIR}/")        
               
rule tigger:
    input:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled.tab',
        IMGT_GERMLINES+'/IMGT_hIGHV.fasta'
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled_distinct-alleles-pass.tab',
        WDIR+'/{patient}_personal_hIGHV.fasta'
    params: name='tigger', partition='long', nproc='4', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{patient}_tigger.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/apply_tigger.R {input_scratch[0]} {input_scratch[1]} {params.nproc} &&\
               cp {output_scratch_concat} {WDIR}/")        
           
rule pool_visits:
    input:
        lambda wildcards: expand(WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE.tab', sample=SAMPLES[wildcards.patient])
    output:
        WDIR+'/{patient}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE_visits-pooled.tab'
    params: name='pool_visits', partition='general', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{patient}_pool_visits.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/pool_visits.R {input_scratch_concat} &&\
               cp {output_scratch} {WDIR}/")         
        
rule filter_functional:
    input:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass.tab'
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-TRUE.tab',
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass_FUNCTIONAL-FALSE.tab'
    params: name='filter_functional', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_filter_functional.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
                set -o nounset &&\
                python {CHANGEO_PATH}/ParseDb.py split -d {input_scratch} -f FUNCTIONAL &&\
                cp {output_scratch_concat} {WDIR}/")
           
rule combine_db_pieces:
    input:
        lambda wildcards: expand(WDIR+'/{piece}_db-pass.tab', piece=IMGT_PIECES[wildcards.sample])
    output:
        WDIR+'/{sample}_minCONSCOUNT2_all-chunks-db-pass.tab'
    params: name='combine_db_pieces', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_combine_db_pieces.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        output_scratch = os.path.basename(output[0])
        shell("cp {input_concat} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               {R_PATH}/bin/Rscript {SCRIPTS}/combine_db_pieces.R {input_scratch_concat} &&\
               cp {output_scratch} {WDIR}/")     
           
rule make_db:
    input:
        WDIR+'/{piece}.fasta',
        WDIR+'/{piece}-extracted'
    output:
        WDIR+'/{piece}_db-pass.tab'
    params: name='make_db', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{piece}_make_db.log'
    run:
        input_scratch = [os.path.basename(s) for s in input]
        output_scratch = os.path.basename(output[0])
        shell("cp {input[0]} $LOCAL_SATA &&\
               cp -r {WDIR}/{wildcards.piece} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate {VIRTUALENV} &&\
               set -o nounset &&\
               python {CHANGEO_PATH}/MakeDb.py imgt -i {wildcards.piece} -s {input_scratch[0]} &&\
               cp {output_scratch} {WDIR}/")
           
rule extract_txz:
    input:
        IMGT_RESULTS+'/{piece}.txz'
    output:
        WDIR+'/{piece}-extracted'
    params: name='extract_txz', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{piece}_extract_txz.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               mkdir {wildcards.piece} &&\
               cd {wildcards.piece} &&\
               tar Jxvf $LOCAL_SATA/{input_scratch} &&\
               cd $LOCAL_SATA &&\
               touch {output_scratch} &&\
               cp -r {wildcards.piece} {WDIR}/ &&\
               cp {output_scratch} {WDIR}/")           
    

    
