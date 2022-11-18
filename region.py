from dataclasses import dataclass, field
from Bio import SeqFeature, Seq
from typing import Optional
import record



@dataclass
class Region:
    gene_name: str
    start: int
    end: int
    strand: int
    record: record.Record
    record_id: int
    feature: SeqFeature.SeqFeature = field(init=False)
    # Needed for the SeedRegion
    start_offset: Optional[int] = field(default=None)
    end_offset: Optional[int] = field(default=None)

    def __post_init__(self):
        self.feature = self.create_feature()

    def create_feature(self) -> SeqFeature.SeqFeature:
        feature_location = SeqFeature.FeatureLocation(self.start, self.end, self.strand)
        return SeqFeature.SeqFeature(feature_location)

    def extract_sequence(self) -> Seq.Seq: 
        return self.feature.extract(self.record.record[self.record_id]).seq
        
    def get_offset(self, start_offset:int, end_offset:int) -> tuple[int, int]:
        genome_len = len(self.record.record[self.record_id])
        if self.strand not in [1, -1]:
            raise ValueError(f"{self.gene_name} lies on unknown strand {self.strand}")
        elif self.strand == 1:
            start = max(0, self.start - start_offset)
            end = min(genome_len, self.start + end_offset)
            return start, end
        elif self.strand == -1:
            start = max(0, self.end - start_offset)
            end = min(genome_len, self.end + end_offset)
            return start, end



class GeneStart(Region):
    # Genomic context for IntaRNA
    FIVE_UTR_OFFSET = 75
    start: int
    end: int

    def __post_init__(self) -> None:
        self.start, self.end = self.get_offset(self.FIVE_UTR_OFFSET, self.FIVE_UTR_OFFSET)
        self.feature = self.create_feature()



class SeedRegion(Region):
    start_offset: int
    end_offset: int

    def __post_init__(self) -> None:
        self.start, self.end = self.get_offset(start_offset=self.start_offset, end_offset=self.end_offset)
        self.feature = self.create_feature()

