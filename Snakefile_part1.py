import os
import glob

##### Paths
VIRTUALENV=             '/local10G/resources/py2.7.9_virtualenv'
CHANGEO_PATH=           '/local10G/debourcy/tools/changeo_0.2.3'
PRESTO_PATH=            '/local10G/debourcy/tools/presto_v0.4.8'
R_PATH=                 '/local10G/resources/R-3.1.2'
IMGT_GERMLINES=         '/local10G/debourcy/tools/tigger_germlines'
SCRIPTS=                '/datastore/debourcy/snakemake_scripts/igh_repertoires'
DATA=                   '/datastore/debourcy/Antibodies/Ellison'
IMGT_RESULTS=           DATA+'/IMGT_output'    
WDIR=                   DATA+'/Snakemake_wdir'
LOGS=                   WDIR+'/logs'
FLU_C_AMPLICONS_REVCOMP= '/local10G/debourcy/Flu_C_amplicons_revcomp.fasta'
R2_PRIMERS=             '/local10G/debourcy/R2_Primers.fasta'
R1_PRIMERS=             '/local10G/debourcy/R1_Primers.fasta'
MUSCLE_PATH=            '/local10G/debourcy/tools/muscle'

if not os.path.exists(WDIR):
    os.makedirs(WDIR)
if not os.path.exists(LOGS):
    os.makedirs(LOGS)

##### Sample lists
SAMPLES = list(set([os.path.basename(FOLDER) for FOLDER in glob.glob(DATA+'/01*')]))

##### Rules
rule all:
    input:
        expand(WDIR+"/{sample}-fastas-for-IMGT-submission",sample=SAMPLES)
        
rule make_chunks:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique_primers-pass_atleast-2.fasta'
    output:
        WDIR+'/{sample}-fastas-for-IMGT-submission'
    params: name='make_chunks', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_make_chunks.log'
    run:
        input_scratch = os.path.basename(input[0])
        output_scratch = os.path.basename(output[0])
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                {R_PATH}/bin/Rscript {SCRIPTS}/make_chunks.R {wildcards.sample} {input_scratch} &&\
                cp *_minCONSCOUNT2*.fasta {WDIR}/ &&\
                touch {output_scratch} &&\
                cp {output_scratch} {WDIR}/")        
        
rule split_conscount:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique_primers-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique_primers-pass_atleast-2.fasta',
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique_primers-pass_under-2.fasta'
    params: name='split_conscount', partition='general', cpus_per_task='4', mem_per_cpu='5300', logfile=LOGS+'/{sample}_split_conscount.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/SplitSeq.py group -s {input_scratch_unzipped} --fasta -f CONSCOUNT --num 2 &&\
                cp {output_scratch_concat} {WDIR}/")

rule confirm_isotypes:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique_primers-pass.fastq.gz',
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique_primers-fail.fastq.gz'
    params: name='confirm_isotypes', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_confirm_isotypes.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/MaskPrimers.py align -s {input_scratch_unzipped} --failed --nproc {params.nproc} --maxerror 0.2 --revpr -p {FLU_C_AMPLICONS_REVCOMP} --maxlen 150 --mode cut --log MP_isotype.out &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")

rule collapse_identicals:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader_collapse-unique.fastq.gz'
    params: name='collapse_identicals', partition='long', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_collapse_identicals.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = os.path.basename(output[0])
        output_scratch_unzipped = os.path.splitext(output_scratch)[0]
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/CollapseSeq.py -s {input_scratch_unzipped} --uf PRCONS --cf CONSCOUNT --act sum --inner &&\
                gzip {output_scratch_unzipped} &&\
                cp {output_scratch} {WDIR}/")

rule consolidate_conscount:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass_reheader.fastq.gz'
    params: name='consolidate_conscount', partition='long', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_consolidate_conscount.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = os.path.basename(output[0])
        output_scratch_unzipped = os.path.splitext(output_scratch)[0]
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/ParseHeaders.py collapse -s {input_scratch_unzipped} -f CONSCOUNT --act min &&\
                gzip {output_scratch_unzipped} &&\
                cp {output_scratch} {WDIR}/")
                
rule filter_missing:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass_missing-pass.fastq.gz'
    params: name='filter_missing', partition='general', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_filter_missing.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/FilterSeq.py missing -s {input_scratch_unzipped} --nproc {params.nproc} -n 10 --inner &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")

rule assemble_pairs:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass.fastq.gz',
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass_assemble-pass.fastq.gz',
        WDIR+'/{sample}_AP.out.gz',
        WDIR+'/{sample}_AP_table.tab.gz'
    params: name='assemble_pairs', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_assemble_pairs.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        input_scratch_unzipped = [os.path.splitext(s)[0] for s in input_scratch]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input_concat} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch_concat} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/AssemblePairs.py align -1 {input_scratch_unzipped[0]} -2 {input_scratch_unzipped[1]} --nproc {params.nproc} --coord presto --rc tail --1f CONSCOUNT --2f CONSCOUNT PRCONS --log {output_scratch_unzipped[1]} --scanrev &&\
                python {PRESTO_PATH}/ParseLog.py -l {output_scratch_unzipped[1]} -f ID OVERLAP ERROR PVALUE &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")
                
rule pair_seq_2:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass.fastq.gz',
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass.fastq.gz',
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass_pair-pass.fastq.gz'
    params: name='pair_seq_2', partition='long', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_pair_seq_2.log'
    run:
        input_concat = ' '.join(input)        
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        input_scratch_unzipped = [os.path.splitext(s)[0] for s in input_scratch]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input_concat} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch_concat} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/PairSeq.py -1 {input_scratch_unzipped[0]} -2 {input_scratch_unzipped[1]} --coord presto &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")                

rule build_R2_consensus:
    input:
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass_align-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass.fastq.gz',
        WDIR+'/{sample}_BC2.out.gz',
        WDIR+'/{sample}_BC2_table.tab.gz'
    params: name='build_R2_consensus', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_build_R2_consensus.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/BuildConsensus.py -s {input_scratch_unzipped} --nproc {params.nproc} --bf BARCODE --pf PRIMER --prcons 0.7 --maxerror 0.1 --log {output_scratch_unzipped[1]} &&\
                python {PRESTO_PATH}/ParseLog.py -l {output_scratch_unzipped[1]} -f BARCODE CONSCOUNT PRIMER DIVERSITY &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")                

rule build_R1_consensus:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass_consensus-pass.fastq.gz',
        WDIR+'/{sample}_BC1.out.gz',
        WDIR+'/{sample}_BC1_table.tab.gz'
    params: name='build_R1_consensus', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_build_R1_consensus.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
				gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/BuildConsensus.py -s {input_scratch_unzipped} --nproc {params.nproc} --bf BARCODE --pf PRIMER --maxerror 0.1 --log {output_scratch_unzipped[1]} &&\
                python {PRESTO_PATH}/ParseLog.py -l {output_scratch_unzipped[1]} -f BARCODE CONSCOUNT PRIMER DIVERSITY &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")                

				
rule align_sets_R2:
    input:
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass_align-pass.fastq.gz'
    params: name='align_sets', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_align_sets_R2.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        r2_primers_scratch = os.path.basename(R2_PRIMERS)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                cp {R2_PRIMERS} $LOCAL_SATA &&\
                python {PRESTO_PATH}/AlignSets.py table -p {r2_primers_scratch} --exec {MUSCLE_PATH} --reverse --outname R2_Primers &&\
                cp R2_Primers_offsets-reverse.tab {WDIR}/ &&\
                python {PRESTO_PATH}/AlignSets.py offset -s {input_scratch_unzipped} --nproc {params.nproc} -d R2_Primers_offsets-reverse.tab &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")                				
				
rule align_sets_R1:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass_align-pass.fastq.gz'
    params: name='align_sets', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_align_sets_R1.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        r1_primers_scratch = os.path.basename(R1_PRIMERS)
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
				cp {R1_PRIMERS} $LOCAL_SATA &&\
                python {PRESTO_PATH}/AlignSets.py table -p {r1_primers_scratch} --exec {MUSCLE_PATH} --reverse --outname R1_Primers &&\
                cp R1_Primers_offsets-reverse.tab {WDIR}/ &&\
                python {PRESTO_PATH}/AlignSets.py offset -s {input_scratch_unzipped} --nproc {params.nproc} -d R1_Primers_offsets-reverse.tab &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")     
               					
rule pair_seq_1:
    input:
        WDIR+'/{sample}_R1_quality-pass_primers-pass.fastq.gz',
        WDIR+'/{sample}_R2_quality-pass_primers-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass_pair-pass.fastq.gz',
        WDIR+'/{sample}_R2_quality-pass_primers-pass_pair-pass.fastq.gz'
    params: name='pair_seq_1', partition='long', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_pair_seq_1.log'
    run:
        input_concat = ' '.join(input)
        input_scratch = [os.path.basename(s) for s in input]
        input_scratch_concat = ' '.join(input_scratch)
        input_scratch_unzipped = [os.path.splitext(s)[0] for s in input_scratch]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        shell("cp {input_concat} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch_concat} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/PairSeq.py -1 {input_scratch_unzipped[0]} -2 {input_scratch_unzipped[1]} --coord illumina --1f BARCODE --2f BARCODE &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")                
                
rule mask_primers_R2:
    input:
        WDIR+'/{sample}_R2_quality-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R2_quality-pass_primers-pass.fastq.gz',
        WDIR+'/{sample}_R2_quality-pass_primers-fail.fastq.gz',
        WDIR+'/{sample}_MP2.out.gz'
    params: name='mask_primers_R2', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_mask_primers_R1.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        UID_LENGTH_R2 = 8
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/MaskPrimers.py score -s {input_scratch_unzipped} --failed --nproc {params.nproc} -p {R2_PRIMERS} --start {UID_LENGTH_R2} --mode cut --barcode --log {output_scratch_unzipped[2]} &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")
                
rule mask_primers_R1:
    input:
        WDIR+'/{sample}_R1_quality-pass.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass_primers-pass.fastq.gz',
        WDIR+'/{sample}_R1_quality-pass_primers-fail.fastq.gz',
        WDIR+'/{sample}_MP1.out.gz'
    params: name='mask_primers_R1', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_mask_primers_R1.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = [os.path.basename(s) for s in output]
        output_scratch_concat = ' '.join(output_scratch)
        output_scratch_unzipped = [os.path.splitext(s)[0] for s in output_scratch]
        output_scratch_unzipped_concat = ' '.join(output_scratch_unzipped)
        UID_LENGTH_R1 = 8
        shell("cp {input} $LOCAL_SATA &&\
                cd $LOCAL_SATA &&\
                gunzip {input_scratch} &&\
                set +o nounset &&\
                source {VIRTUALENV}/bin/activate &&\
                set -o nounset &&\
                python {PRESTO_PATH}/MaskPrimers.py score -s {input_scratch_unzipped} --failed --nproc {params.nproc} -p {R1_PRIMERS} --start {UID_LENGTH_R1} --mode cut --barcode --log {output_scratch_unzipped[2]} &&\
                gzip {output_scratch_unzipped_concat} &&\
                cp {output_scratch_concat} {WDIR}/")


rule filter_quality_R2:
    input:
        DATA+'/{sample}/{sample}_R2.fastq.gz'
    output:
        WDIR+'/{sample}_R2_quality-pass.fastq.gz'
    params: name='filter_quality_R2', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_filter_quality_R2.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = os.path.basename(output[0])
        output_scratch_unzipped = os.path.splitext(output_scratch)[0]
        shell("cp {input} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               gunzip {input_scratch} &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate &&\
               set -o nounset &&\
               python {PRESTO_PATH}/FilterSeq.py quality -s {input_scratch_unzipped} --nproc {params.nproc} -q 20 &&\
               gzip {output_scratch_unzipped} &&\
               cp {output_scratch} {WDIR}/")

rule filter_quality_R1:
    input:
        DATA+'/{sample}/{sample}_R1.fastq.gz'
    output:
        WDIR+'/{sample}_R1_quality-pass.fastq.gz'
    params: name='filter_quality_R1', partition='long', nproc='12', cpus_per_task='12', mem_per_cpu='5300', logfile=LOGS+'/{sample}_filter_quality_R1.log'
    run:
        input_scratch = os.path.basename(input[0])
        input_scratch_unzipped = os.path.splitext(input_scratch)[0]
        output_scratch = os.path.basename(output[0])
        output_scratch_unzipped = os.path.splitext(output_scratch)[0]
        shell("cp {input} $LOCAL_SATA &&\
               cd $LOCAL_SATA &&\
               gunzip {input_scratch} &&\
               set +o nounset &&\
               source {VIRTUALENV}/bin/activate &&\
               set -o nounset &&\
               python {PRESTO_PATH}/FilterSeq.py quality -s {input_scratch_unzipped} --nproc {params.nproc} -q 20 &&\
               gzip {output_scratch_unzipped} &&\
               cp {output_scratch} {WDIR}/")


    
