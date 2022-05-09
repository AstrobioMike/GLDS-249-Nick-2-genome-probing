####################################################################
##### Snakefile for Lee work on GLDS-249 Nick 2-genome probing #####
####################################################################

import os

configfile: "config.yaml"


########################################
############# General Info #############
########################################

"""
See the corresponding 'config.yaml' file for general use information.
Variables that may need to be adjusted should be changed there.
"""

########################################
#### Reading samples file into list ####
########################################

genome_IDs_file = config["genome_IDs_file"]
genome_ID_list = [line.strip() for line in open(genome_IDs_file)]

sample_IDs_file = config["sample_IDs_file"]
sample_ID_list = [line.strip() for line in open(sample_IDs_file)]

########################################
######## Setting up directories ########
########################################


dirs_to_create = [
                  config["input_anvio_files_dir"], 
                  config["bam_files_dir"], 
                  config["bowtie2_indexes_dir"], 
                  config["contigs_dbs_dir"], 
                  config["profile_dbs_dir"],
                  config["logs_dir"]
                  ]

for directory in dirs_to_create:
    try:
        os.mkdir(directory)
    except:
        pass


# making lists of profiles that will need to be merged
profiles_dict = {}
for genome in genome_ID_list:

    curr_list_of_profiles = expand(config["profile_dbs_dir"] + "{sample}-to-" + str(genome) + "-profile/PROFILE.db", sample = sample_ID_list)

    profiles_dict[genome] = curr_list_of_profiles

########################################
############# Rules start ##############
########################################


rule all:
    input:
        anvio_input_files = expand(os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa"), genome = genome_ID_list),
        contigs_dbs = expand(os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db"), genome = genome_ID_list),
        bowtie2_indexes = expand(os.path.join(config["bowtie2_indexes_dir"], "{genome}.1.bt2"), genome = genome_ID_list),
        bam_files = expand(os.path.join(config["bam_files_dir"] + "{sample}-to-{genome}.bam"), sample = sample_ID_list, genome = genome_ID_list),
        profile_dbs = expand(os.path.join(config["profile_dbs_dir"], "{sample}-to-{genome}-profile/PROFILE.db"), sample = sample_ID_list, genome = genome_ID_list),
        merged_profiles = expand(os.path.join(config["profile_dbs_dir"], "{genome}-merged-profile"), genome = genome_ID_list)


rule genbank_to_anvio:
    input:
        genbank = os.path.join(config["genbank_files_dir"], "{genome}.gb")
    output:
        fasta = os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa"),
        gene_calls = os.path.join(config["input_anvio_files_dir"], "{genome}-external-gene-calls.txt"),
        gene_functions = os.path.join(config["input_anvio_files_dir"], "{genome}-external-functions.txt")
    params:
        output_prefix = os.path.join(config["input_anvio_files_dir"], "{genome}")
    log:
        config["logs_dir"] + "genbank_to_anvio-{genome}.log"
    shell:
        """
        anvi-script-process-genbank -i {input.genbank} -O {params.output_prefix} > {log} 2>&1
        """


rule make_and_annotate_contigs_db:
    input:
        fasta = os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa"),
        gene_calls = os.path.join(config["input_anvio_files_dir"], "{genome}-external-gene-calls.txt"),
        gene_functions = os.path.join(config["input_anvio_files_dir"], "{genome}-external-functions.txt"),
    output:
        contigs_db = os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db")
    params:
        name = "{genome}",
        num_threads = config["general_anvio_threads"],
        COGs_dir = config["anvio_COGs_dir"],
        KOs_dir = config["anvio_KOs_dir"]
    log:
        config["logs_dir"] + "make_and_annotate_contigs_db-{genome}.log"
    shell:
        """
        anvi-gen-contigs-database -f {input.fasta} -o {output.contigs_db} -n {params.name} --external-gene-calls {input.gene_calls} \
                                  -T {params.num_threads} --ignore-internal-stop-codons > {log} 2>&1
                                  
        anvi-import-functions -c {output.contigs_db} -i {input.gene_functions} >> {log} 2>&1

        anvi-run-hmms -T {params.num_threads} -I Bacteria_71 -c {output.contigs_db} >> {log} 2>&1

        anvi-scan-trnas -T {params.num_threads} -c {output.contigs_db} >> {log} 2>&1

        anvi-run-ncbi-cogs -c {output.contigs_db} --cog-data-dir {params.COGs_dir} -T {params.num_threads} --sensitive >> {log} 2>&1

        anvi-run-kegg-kofams -c {output.contigs_db} --kegg-data-dir {params.KOs_dir} -T {params.num_threads} >> {log} 2>&1
        """


rule make_bowtie2_index:
    input:
        fasta = os.path.join(config["input_anvio_files_dir"], "{genome}-contigs.fa")
    output:
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.1.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.2.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.3.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.4.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.rev.1.bt2"),
        os.path.join(config["bowtie2_indexes_dir"], "{genome}.rev.2.bt2")
    params:
        output_prefix = os.path.join(config["bowtie2_indexes_dir"], "{genome}"),
        num_threads = config["mapping_threads"]
    log:
        config["logs_dir"] + "make_bowtie2_index-{genome}.log"
    shell:
        """
        bowtie2-build --threads {params.num_threads} {input.fasta} {params.output_prefix} > {log} 2>&1
        """


rule mapping:
    input:
        index_built_trigger = os.path.join(config["bowtie2_indexes_dir"], "{genome}.1.bt2"),
        R1 = os.path.join(config["reads_dir"] + "{sample}" + config["read_R1_suffix"]),
        R2 = os.path.join(config["reads_dir"] + "{sample}" + config["read_R1_suffix"])
    output:
        bam = os.path.join(config["bam_files_dir"] + "{sample}-to-{genome}.bam")
    params:
        index = os.path.join(config["bowtie2_indexes_dir"], "{genome}"),
        num_threads = config["mapping_threads"]
    log:
        os.path.join(config["bam_files_dir"] + "{sample}-to-{genome}-mapping-info.txt")
    shell:
        """
        bowtie2 -q --threads {params.num_threads} -x {params.index} -1 {input.R1} -2 {input.R2} --no-unal 2> {log} | samtools view -b | samtools sort -@ {params.num_threads} > {output.bam}
        samtools index -@ {params.num_threads} {output.bam}
        """


rule anvi_profile:
    input:
        contigs_db = os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db"),
        bam = os.path.join(config["bam_files_dir"] + "{sample}-to-{genome}.bam")
    output:
        profile_db = os.path.join(config["profile_dbs_dir"], "{sample}-to-{genome}-profile/PROFILE.db")
    params:
        profile_db_dir = os.path.join(config["profile_dbs_dir"], "{sample}-to-{genome}-profile"),
        num_threads = config["general_anvio_threads"]
    log:
        log = config["logs_dir"] + "anvi_profile-{sample}-to-{genome}.log"
    shell:
        """
        name=$(printf {wildcards.sample} | tr "[\-.]" "_")
        anvi-profile -c {input.contigs_db} -i {input.bam} -o {params.profile_db_dir} -S ${{name}} \
                     --cluster-contigs --min-contig-length 500 -T {params.num_threads} \
                     --overwrite-output-destinations > {log} 2>&1
        """


rule anvi_merge:
    input:
        contigs_db = os.path.join(config["contigs_dbs_dir"], "{genome}-contigs.db"),
        profiles = lambda wildcards: profiles_dict[wildcards.genome]
    output:
        merged_profile = os.path.join(config["profile_dbs_dir"], "{genome}-merged-profile")
    log:
        log = config["logs_dir"] + "anvi_merge-{genome}.log"
    shell:
        """
        name=$(printf {wildcards.genome} | tr "[\-.]" "_")

        anvi-merge -c {input.contigs_db} -o {output.merged_profile} -S ${{name}} --skip-hierarchical-clustering {input.profiles} > {log} 2>&1
        
        # getting split default order
        anvi-export-table --table splits_basic_info -f split {input.contigs_db} -o {wildcards.genome}-split-order.tmp

        tail -n +2 {wildcards.genome}-split-order.tmp > {wildcards.genome}-split-order.txt
        rm {wildcards.genome}-split-order.tmp

        # adding to merged profile db
        anvi-import-items-order -i {wildcards.genome}-split-order.txt -p {output.merged_profile}/PROFILE.db --make-default --name default
        """


rule clean_all:
    shell:
        """
        rm -rf {dirs_to_create}
        """
