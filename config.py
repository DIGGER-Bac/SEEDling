import yaml
from dataclasses import dataclass, field


@dataclass
class Config():
    config_path: str 
    config: dict = field(init=False)

    subject_path: str = field(init=False)
    target_path: str = field(init=False)
    output_path: str = field(init=False)
    select_top: int = field(init=False)
    start_offset: int = field(init=False)
    end_offset: int = field(init=False)
    step_size: int = field(init=False)
    seq_length: int = field(init=False)
    srna_prefix: str = field(init=False, default="")
    srna_suffix: str = field(init=False, default="")
    srna_template: str = field(init=False)
    exclude_sequences_path: str = field(init=False)
    include_genes_path: str = field(init=False)
    blast_evalue: float = field(init=False)
    
    def __post_init__(self) -> None:
        with open(self.config_path) as file:
            self.config = yaml.load(file, Loader=yaml.FullLoader)
        self.load_parameters()

    def load_parameters(self) -> None:
        self.subject_path = self.config["subject_path"]
        self.target_path = self.config["target_path"]
        self.output_path = self.config["output_path"]
        self.select_top = self.config["select_top"]
        self.start_offset = self.config["start_offset"]
        self.end_offset = self.config["end_offset"]
        self.step_size = self.config["step_size"]
        self.seq_length = self.config["seq_length"]
        self.srna_prefix = self.config["srna_prefix"]
        self.srna_suffix = self.config["srna_suffix"]
        self.srna_template = self.config["srna_template"]
        self.exclude_sequences_path = self.config["exclude_sequences_path"]
        self.include_genes_path = self.config["include_genes_path"]
        self.blast_evalue = self.config["blast_evalue"]