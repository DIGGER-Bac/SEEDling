# SEEDling 
## Description
SEEDling predicts the best possible synthetic sRNA sequence for the target(s) of choice. In order to do so it requires the following user inputs: 
1.  an annotated target genbank file (single gene to
whole genome)
2.  a wild-type sRNA and its corresponding scaffold
3.  the reference genome of the target organism and 
4.  a config file    

Based on the input files SEEDling predicts the best antisense regions, uses [BLASTn](https://pubmed.ncbi.nlm.nih.gov/20003500/) to determine off-target effects by the expect value (E-value). Afterwards, it combines the seed region with the provided sRNA scaffold to
determine the binding energy via [intaRNA](https://pubmed.ncbi.nlm.nih.gov/28472523/) and finally applies [RNApdist](https://pubmed.ncbi.nlm.nih.gov/22115189/) to compare synthetic sRNAs and wild-type sRNA in regard of potential folding alterations, allowing to exclude synthetic sRNAs with drastic structural changes.

## Requirements
The tool was tested running Ubuntu 22.04.1 LTS with the following installations:   
- [Python 3.9.13](https://www.python.org/downloads/release/python-3913/)   
- [conda 4.13.0](https://github.com/conda/conda/releases/tag/4.13.0)    
- [IntaRNA 3.3.1 ](https://github.com/BackofenLab/IntaRNA/releases/tag/v3.3.1)    
- [ViennaRNA 2.4.17](https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.4.17)  
- [NCBI-Blast 2.12.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/)
- [PyYAML 6.0](https://pyyaml.org/)
- [Biopython 1.79](https://biopython.org/wiki/Download)
- [pandas 1.5.1](https://pandas.pydata.org/)



## Usage
## Docker
Due to the number of dependencies, it is recommended to run the pipeline via the Docker image we've created. While docker can be run from the command line interface, we recommend installing the [Docker Desktop](https://www.docker.com/products/docker-desktop/) app.
If you're new to Docker follow the quick start guide that can be found [here](https://docs.docker.com/desktop/get-started/).
The SEEDling Docker image can be downloaded from [here](https://owncloud.gwdg.de/index.php/s/7tKXsNXfq9OdQzs). To load the image please follow the instructions from the official [Docker documentation](https://docs.docker.com/engine/reference/commandline/load/).

### Creating the Docker Container
When creating the container, it is important to set the following optional settings:   

**Volumes:**
> **Host path**: File path to an empty folder on your host machine.

> **Container path**: /home/DIGGER/SEEDling/input

Providing these settings will allow you to interact with the container by parsing files into the host folder. 

### Running SEEDling in Docker
1. Start and connect to the container.
2. Switch the shell from `sh` to `bash` by entering and executing the `bash` command.
3. Navigate to the script directory by executing 
```cd home/DIGGER/SEEDling```
4. Download the example *input* files from GitHub and place the files into the *host path* directory.
5. Change the files as needed (Note: Make sure to change all paths in the config.yml to include "input/" as a prefix)
6. Execute SEEDling by running
```python3 SEEDling.py -c input/config.yml```
7. SEEDling now executes and generates an output file as defined in the *config.yml*

## Running SEEDling without Docker
Before running SEEDling, install all dependencies that are listed in the *requirement* section. Make sure to add *ViennaRNA*, *IntaRNA* and *BLAST* to your *PATH*-variable. Change the configuration file to the the parameters of your choice. 

After setting everything up, SEEDling can be run as follows:
> python3 SEEDling.py -c input/config.yml

## Configuration File
The following parameters must be set in the configuration [YAML](https://yaml.org/) file:   
(If not needed, you can leave the parameter empty)
| Parameter | Type | Example | Description |
|-----------|------|---------|-------------|
| subject_path   | string | "input/e_coli_k12_mg1655.gb" | Path to GenBank target genome |
| target_path | string | "input/acrAB_operon.gbk" |Path to GenBank file of genes for which the SEED will be calculated |
| output_path | string | output.csv | Path to output file (.csv) |
| select_top | integer | 10 | Number of SEED predictions per gene |
| start_offset | integer | 0 | Offset (left of the start codon)|
| end_offset | integer | 0 | Offset (right of the start codon)|
| stepsize | integer | 1 | Step size (for the sliding window)
| seq_length | integer | 16 | Length of antisense sequence |
| seq_prefix | string | GCTA | Prefix (scaffold) |
| seq_suffix | string | ATCG | Suffix (scaffold) |
| srna_template | string | GCCACTGCTTTTCTTTGATGTCCCCATTTTG | Template scaffold |
| exclude_sequences_path | string | input/exclude_sites.fasta | Path to FASTA file that contains sequences to exclude (i.e. restriction recognition sites) |
| include_genes_path | string | input/filter_short.txt | Path to *.txt file of gene names that should be chosen (If left empty all genes will be chosen)| 
| blast_evalue | float | 0.04 | E-value for BLASTn offtarget checks |

### Output
Following output file be written to the specified output *.csv file.    

| Parameter | Example Output | Description |
|-----------|----------------|-------------|
| ID | ULR0V2S0 | Unique sRNA identifier |
| Source | input/recA_operon.gbk | Input file path |
| RNAdist | 10.8959 | Distance to the template structure | 
| Hybrid | -41.69 | RNA hybridization energy (kcal/mol) | 
| Prefix | AACAG | Prefix sequence as defined in the input |
| Suffix | ATCAC | Suffix sequence as defined in the input |
| Seed | CCTGGCTCATCATACGTGCCGC | Antisense sequence generated by to tool |
| Gene | recA | Gene name of the gene corresponding to the current antisense |
| Start | 580 | Starting position (in bp) of the antisense sequence |
| End | 602 | End position (in bp) of the antisense sequence|
| Strand | + | Strand of the current gene |
| Offsite | False | Found offtargets (boolean) |
| Fold | ..(((.....))).... | Secondary structure of the scaffold + seed |
| FullSeq | AACAGCCTGGCTCATCATACGTGCCGCATCAC | Prefix + seed + suffix sequence |



## References
Martin Mann, Patrick R. Wright, Rolf Backofen, IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions, Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W435–W439, https://doi.org/10.1093/nar/gkx279   

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. In BMC Bioinformatics (Vol. 10, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1471-2105-10-421

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, 1 June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163

Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. In Algorithms for Molecular Biology (Vol. 6, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1748-7188-6-26