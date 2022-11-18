import os
import pandas as pd 
import logging 
import record
from config import Config
from gene_target import GeneTarget
from srna import sRNA
from blast import Blast
from region import SeedRegion, GeneStart



def main(config_path:str) -> None:
    logging.info("Loading configurations ...")
    config: Config = Config(config_path=config_path)
    target_record: record.Target = record.Target(config.target_path)
    logging.info(f"Successfully loaded target {config.target_path}")
    subject_record: record.Subject = record.Subject(config.subject_path)            
    logging.info(f"Successfully loaded subject {config.subject_path}")


    # Get target gene names
    with open(config.include_genes_path) as file:
        read_buffer = file.read()
    target_gene_names:list[str] = [gn for gn in read_buffer.split("\n") if gn]

    # Get target genes
    target_genes: list[GeneTarget] = [GeneTarget(gn, config, target_record) for gn in target_gene_names]

    seed_regions: list[SeedRegion] = []
    genomic_context:dict[str, GeneStart] = {}
    for target_gene in target_genes:
        logging.info(f"Processing {target_gene.gene_name}")
        genomic_context[target_gene.gene_name] = GeneStart(
            gene_name=target_gene.gene_name,
            start=target_gene.gene_feature.location.start,
            end=target_gene.gene_feature.location.end,
            strand=target_gene.gene_feature.location.strand,
            record=target_gene.target_record,
            record_id=target_gene.record_id
        )

        for i in range(config.start_offset + config.end_offset - config.seq_length):
            seed_start = config.start_offset - i
            seed_end = -seed_start+config.seq_length
            seed_regions.append(SeedRegion(
                gene_name = target_gene.gene_name,
                start = target_gene.gene_feature.location.start,
                end = target_gene.gene_feature.location.end,
                strand = target_gene.gene_feature.location.strand,
                record = target_record,
                record_id = target_gene.record_id,
                start_offset = seed_start,
                end_offset = seed_end
            ))


    srnas:list[sRNA] = []
    for seed in seed_regions:
        srnas.append(sRNA(
            seed = seed,
            context = genomic_context[seed.gene_name],
            config = config
        ))

    # Run BLASTn
    generate_blast_result(config, srnas, subject_record)
    # Perform selection
    final_selection: pd.DataFrame =  perform_selection(srnas, config.select_top)
    # Write results to file
    write_output(final_selection, config.output_path)

    # Cleanup
    try:
        os.remove("1_dp.ps")
        os.remove("2_dp.ps")
    except FileNotFoundError as e:
        logging.warning("Unable to cleanup 1/2_dp.ps")
        pass


def generate_blast_result(config:Config, srnas:list[sRNA], subject_record:record.Subject) -> None:
    blast = Blast(config=config, srnas=srnas, subject=subject_record)
    for srna in srnas:
        if srna.identifier in blast.offsite_found:
            srna.blast_offsite = True


def perform_selection(srnas:list[sRNA], select_top:int) -> pd.DataFrame:
    data = [s.get_values() for s in srnas]
    
    header = ['ID', 'Source', 'RNAdist', 'HybridEnergy', 'Prefix', 'Suffix', 'Seed', 'Gene', 'Start', 'End', 'Strand', 'Offsite', 'IllegalSite','Fold', 'FullSeq']
    df = pd.DataFrame(data=data, columns=header)
    grouped = df.groupby("Gene")

    results_df = []
    for _, group in grouped:
        filtered = group[(group["Offsite"]==False) & (group["IllegalSite"]==False)]
        filtered = filtered.sort_values(by=["RNAdist"])
        filtered["RNA_score"] = list(range(len(filtered)))
        filtered = filtered.sort_values(by=["HybridEnergy"]) 
        filtered["Energy_score"] = list(range(len(filtered)))
        filtered["Total_score"] = filtered["RNA_score"] + filtered['Energy_score']
        results_df.append(filtered.sort_values(by=["Total_score"])[:select_top])
        
    return pd.concat(results_df)


def write_output(final_selection: pd.DataFrame, output_path:str) -> None:
    logging.info("Writing output...")
    final_df = final_selection.drop(['RNA_score', 'Energy_score', 'Total_score'], axis=1)
    final_df.to_csv(output_path, index=False)
    