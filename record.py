import os
from dataclasses import dataclass, field
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO



@dataclass
class Record:
    path:str
    file_type: str = field(init=False)
    record: list[SeqRecord] = field(init=False)

    def __post_init__(self):
        self.file_type = self.get_file_type()
        self.record = self.get_record()


    def get_record(self) -> list[SeqRecord]:
        records = []
        for record in SeqIO.parse(self.path, self.file_type):
            records.append(record)

        return records
    

    def get_file_type(self) -> str:
        filename = os.path.basename(self.path)
        fileformat = filename.split(".")[-1]
        if fileformat in ["mfa", "fa", "fasta"]:
            return "fasta"
        elif fileformat in ["gb", "gbk"]:
            return "genbank"
        else:
            raise ValueError(f"Unknown file extension {fileformat}")


@dataclass
class Target(Record):
    pass


@dataclass
class Subject(Record):
    pass


@dataclass
class Exclude(Record):
    pass
