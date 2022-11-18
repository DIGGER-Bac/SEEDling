from dataclasses import dataclass, field
from Bio.SeqFeature import SeqFeature
import record
from config import Config
import logging



@dataclass
class GeneTarget:
    gene_name: str
    config: Config
    target_record: record.Target

    gene_feature: SeqFeature = field(init=False)
    record_id: int = field(init=False)
    
    def __post_init__(self) -> None:
        self.record_id, self.gene_feature = self.select_feature()
        logging.info(f"Gene {self.gene_name} found")
        pass

    def select_feature(self) -> tuple[int, SeqFeature]:
        for rec_id, record in enumerate(self.target_record.record):
            for feat in record.features:
                if feat.type != "CDS":
                    continue
                try:
                    if feat.qualifiers["gene"][0] == self.gene_name:
                        return rec_id, feat
                except KeyError as e:
                    pass

                try:
                    if feat.qualifiers["label"][0] == self.gene_name:
                        return rec_id, feat
                except KeyError as e:
                    pass
        
        error_msg = f"Unable to find the gene {self.gene_name}"
        logging.error(f"Unable to find the gene {self.gene_name}")
        raise ValueError(error_msg)