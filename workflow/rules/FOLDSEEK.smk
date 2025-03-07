configfile: "config/config.yaml"

import os
include: "AF3.smk"

OUTPUT_DIR = config["output_dir"] # check if can be imported from AF3.smk
FOLDSEEK_CONTAINER = config["foldseek_flags"]["--foldseek_container"]
#TOP_N_HITS = config["top_n_hits"]
FOLDSEEK_DBS_DIR = config["foldseek_flags"]["--foldseek_databases_dir"]
#os.makedirs(FOLDSEEK_DBS_DIR,exist_ok=True)

def get_foldseek_outputs(wildcards):
    PREPROCESSING_DIR = checkpoints.PREPROCESSING.get(**wildcards).output[0]
    JOB_NAMES, = glob_wildcards(os.path.join(PREPROCESSING_DIR, "{i}.json"))
    return  list(expand(os.path.join(OUTPUT_DIR,"FOLDSEEK_EASY_SEARCH","{i}_against_{db}.tsv"),i=JOB_NAMES,db=["alphafold_uniprot",
                                                                                                        "alphafold_uniprot50-minimal",
                                                                                                        "alphafold_uniprot50",
                                                                                                        "alphafold_proteome",
                                                                                                        "alphafold_swiss-prot",
                                                                                                        "esmatlas30",
                                                                                                        "pdb",
                                                                                                        "cath50",
                                                                                                        "bfmd",
                                                                                                        "bfvd",
                                                                                                        "prostt5",
                                                                                                        ]))

rule foldseek_all:
    input:
        get_foldseek_outputs,
        os.path.join(OUTPUT_DIR,"FOLDSEEK_DOWNLOAD_DBS","download_db_done.txt")

rule FOLDSEEK_DOWNLOAD_DBS:
    container:
        FOLDSEEK_CONTAINER
    params:
        db_dir = FOLDSEEK_DBS_DIR,
    output:
        touch(os.path.join(OUTPUT_DIR,"FOLDSEEK_DOWNLOAD_DBS","download_db_done.txt"))
    shell:
        """
        mkdir -p output/FOLDSEEK_DOWNLOAD_DBS/dbs
        foldseek_avx2 databases BFVD output/FOLDSEEK_DOWNLOAD_DBS/dbs/bfvd tmpDir
        """

# foldseek_avx2 databases Alphafold/UniProt {params.db_dir}/alphafold_uniprot tmpDir
# foldseek_avx2 databases Alphafold/UniProt50-minimal {params.db_dir}/alphafold_uniprot50-minimal tmpDir
# foldseek_avx2 databases Alphafold/UniProt50 {params.db_dir}/alphafold_uniprot50 tmpDir
# foldseek_avx2 databases Alphafold/Proteome {params.db_dir}/alphafold_proteome tmpDir,
# foldseek_avx2 databases Alphafold/Swiss-Prot {params.db_dir}/alphafold_swiss-prot tmpDir
# foldseek_avx2 databases ESMAtlas30 {params.db_dir}/esmatlas30 tmpDir
# foldseek_avx2 databases PDB {params.db_dir}/pdb tmpDir
# foldseek_avx2 databases CATH50 {params.db_dir}/cath50 tmpDir
# foldseek_avx2 databases BFMD {params.db_dir}/bfmd tmpDir
# foldseek_avx2 databases ProstT5 {params.db_dir}/prostt5 tmpDir
rule FOLDSEEK_EASY_SEARCH_AGAINST_ALPHAFOLD_UNIPROT:
    input:
        os.path.join(OUTPUT_DIR,"AF3_INFERENCE","{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"alphafold_uniprot")
    output:
        os.path.join(OUTPUT_DIR,"FOLDSEEK_EASY_SEARCH","{i}_against_alphafold_uniprot.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """


rule FOLDSEEK_EASY_SEARCH_AGAINST_ALPHAFOLD_UNIPROT50_MINIMAL:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"alphafold_uniprot50-minimal")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_alphafold_uniprot50-minimal.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_ALPHAFOLD_UNIPROT50:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"alphafold_uniprot50")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_alphafold_uniprot50.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_ALPHAFOLD_PROTEOME:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"alphafold_proteome")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_alphafold_proteome.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_ALPHAFOLD_SWISS_PROT:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"alphafold_swiss-prot")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_alphafold_swiss-prot.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_ESMATLAS30:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"esmatlas30")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_esmatlas30.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_PDB:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"pdb")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_pdb.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_CATH50:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"cath50")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_cath50.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_BFMD:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"bfmd")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_bfmd.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_BFVD:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"bfvd")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_bfvd.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """

rule FOLDSEEK_EASY_SEARCH_AGAINST_PROSTT5:
    input:
        os.path.join(OUTPUT_DIR, "AF3_INFERENCE", "{i}/{i}/{i}_model.cif"),
    params:
        db = os.path.join(FOLDSEEK_DBS_DIR,"prostt5")
    output:
        os.path.join(OUTPUT_DIR, "FOLDSEEK_EASY_SEARCH", "{i}_against_prostt5.tsv"),
    container:
        FOLDSEEK_CONTAINER
    shell:
        """
        foldseek_avx2 \
        easy-search \
        {input} \
        {params.db} {output} tmp \
        --alignment-type 1 \
        --format-output "query,target,complexqtmscore,complexttmscore,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,complexassignid"
        """
