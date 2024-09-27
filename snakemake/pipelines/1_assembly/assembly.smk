import pandas as pd

configfile: "config/_______.yaml"

samples_df = pd.read_csv("tsv/______.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        

rule fastqc:
    """
    
    """
    input:
        
    output:
        
    conda:
        
    params:
       
    shell:
        """
       
        """

rule trim_galore:
    """
    
    """
    input:
        
    output:
        
    conda:
        
    params:
       
    shell:
        """
       
        """