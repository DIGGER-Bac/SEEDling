import re
from subprocess import Popen, PIPE
from dataclasses import dataclass, field
import random
import string
from Bio import Seq

from record import Exclude
from config import Config
from region import SeedRegion, GeneStart



@dataclass
class sRNA:
    seed: SeedRegion
    context: GeneStart
    config: Config

    # reverse_complement() of the region sequence
    sequence: Seq.Seq = field(init=False) 
    scaffold: Seq.Seq = field(init=False)
    secondary_structure: str = field(init=False)
    contains_illegal_site: bool = field(init=False)
    secondary_struct_distance: float = field(init=False)
    identifier: str = field(init=False)
    hybrid_energy: float = field(init=False)
    blast_offsite: bool = field(init=False, default=False)

    def __post_init__(self) -> None:
        sequence_buffer = self.seed.extract_sequence()
        self.sequence = sequence_buffer.reverse_complement()
        self.scaffold =  self.get_scaffold()
        self.secondary_structure = self.calc_sec_structure()
        self.contains_illegal_site = self.check_illegal_site()
        self.secondary_struct_distance = self.calc_rna_pdist(self.config.srna_template)
        self.identifier = ''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
        self.hybrid_energy = self.calc_hybrid_energy()


    def get_scaffold(self) -> Seq.Seq:
        # appends the prefix and suffix to the sRNA seed
        prefix = self.config.srna_prefix
        suffix = self.config.srna_suffix
        return Seq.Seq(prefix) + self.sequence + Seq.Seq(suffix)


    def check_illegal_site(self) -> bool:
        exclude = Exclude(self.config.exclude_sequences_path)
        for rec in exclude.record:
            if (rec.seq in self.scaffold) or (rec.reverse_complement().seq in self.scaffold):
                return True
        return False


    def calc_hybrid_energy(self) -> float:
        query = str(self.context.extract_sequence())
        cmd = f"IntaRNA -t {self.scaffold} -q {query} --outmode=C --outCsvCols=start1,end1,start2,end2,E_hybrid"
        result = self.run_cmd(cmd)

        csv_data = result.split("\n")
        headers = ["start1", "end1", "start2", "end2", "E_hybrid"]

        try:
            headers = csv_data[0].split(";")
            value_list = csv_data[1].split(";")
        except Exception as e:
            print(csv_data)
            print(result)
            print(e)
            exit()

        try:
            return_dict = dict(zip(headers, value_list))
            ret_value = return_dict["E_hybrid"]
            return float(ret_value)
        except KeyError as e:
            return 0.0


    def calc_sec_structure(self) -> str:
        # Calculates the secondary structure of the full scaffold sRNA
        cmd = f"echo {str(self.scaffold)} | RNAfold --noPS /dev/stdin"
        result = self.run_cmd(cmd)
        pattern = r"[().]+"
        struct = re.findall(pattern, result)[0]

        return str(struct)


    def __str__(self) -> str:
        return ",".join([str(osl) for osl in self.get_values()])


    def get_values(self) -> list:
        out_list = [
            self.identifier,    # Sequence identifier
            self.config.target_path,    # Path to target file
            self.secondary_struct_distance, # Secondary structure distance to template
            float(self.hybrid_energy), # IntaRNA hybrid energy value
            str(self.config.srna_prefix),    # defined srna prefix
            str(self.config.srna_suffix),    # defined srna suffix
            self.sequence,   # seed region
            self.seed.gene_name,  # target gene name
            self.seed.start,    # start position of seed region
            self.seed.end,      # end position of seed region
            self.seed.strand,   # strand of seed region
            self.blast_offsite, # bool if offsites where found
            self.contains_illegal_site,
            self.secondary_structure,   # rna secondary structure of scaffold (vienna)
            str(self.scaffold)   # seed + pre/suffix sequences
        ]

        return out_list


    def run_cmd(self, cmd: str, inp: str = "") -> str:
        process = Popen(args=cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = process.communicate(input=inp.encode())
        if stderr:
            print(stderr.decode("ascii"))
        return stdout.decode("ascii")


    def calc_rna_pdist(self, reference:str) -> float:
        if not reference:
            return 0.0
        process = Popen(args="RNApdist", stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = process.communicate(input=f"{self.scaffold}\n{reference}\n".encode())

        if stderr:
            print(stderr.decode("ascii"))

        return float(stdout.decode("ascii"))