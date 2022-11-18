import os
import random
import string
from dataclasses import dataclass, field
from subprocess import Popen, PIPE
from collections import Counter

from Bio import SeqIO

from config import Config
from srna import sRNA
from record import Subject




@dataclass
class Blast:
    config: Config
    srnas: list[sRNA]
    subject: Subject

    query_filename: str = field(init=False)
    subject_filename: str = field(init=False)
    result_filename: str = field(init=False)

    blast_result: list[dict] = field(init=False) 
    offsite_found: list[str] = field(init=False)

    def __post_init__(self) -> None:
        self.query_filename = self.generate_name()
        self.subject_filename = self.generate_name()
        self.result_filename = self.generate_name()
        self.write_query()
        self.run_blast()
        self.blast_result = self.parse_results()
        self.cleanup()
        self.offsite_found = self.filter_results()


    def generate_name(self) -> str:
        return ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))


    def write_query(self) -> None:
        with open(self.subject_filename, "w") as file:
            SeqIO.write(self.subject.record, file, "fasta")

        with open(self.query_filename, "w") as file:
            for srna in self.srnas:
                file.write(f">{srna.identifier}\n") 
                file.write(f"{str(srna.scaffold)}\n") 


    def run_blast(self) -> None:
        cmd = f"blastn -query {self.query_filename} -subject {self.subject_filename} -outfmt 6 -out {self.result_filename} -evalue {self.config.blast_evalue}" # -word_size {int(config.seq_length/2)}"
        process = Popen(args=cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
        _, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))


    def parse_results(self) -> list[dict]:
        with open(self.result_filename) as file:
            data = file.read()

        headers = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
        lines = [ret for ret in data.split("\n") if ret]
        hits = [b.split("\t") for b in lines]

        result = []
        for br in hits:
            br_dict = dict(zip(headers, br))
            result.append(br_dict)

        return result

    
    def filter_results(self) -> list[str]:
        qseqids =  [bl["qseqid"] for bl in self.blast_result]
        qseq_counter = Counter(qseqids)
        return [x[0] for x in qseq_counter.items() if x[1]>=2]


    def cleanup(self) -> None:
        try:
            os.remove(self.subject_filename)
            os.remove(self.query_filename)
            os.remove(self.result_filename)
        except FileNotFoundError as e: 
            # do some logging
            pass