# SEEDling 
## Description
SEEDling predicts the best possible synthetic sRNA sequence for the desired target(s). In order to do this, it requires the following user input:
1.  Annotated target GenBank file (single gene to
whole genome)
2.  Wild-type sRNA and its corresponding scaffold
3.  Reference genome of the target organism and 
4.  Config file (see below)    


Based on these inputs, SEEDling iterates over the range of target regions to find the best candidate SEED regions. This is done by combining the specified sRNA prefix/suffix sequences with the respective SEED sequence to form an sRNA scaffold of which the binding energy is determined using [intaRNA](https://pubmed.ncbi.nlm.nih.gov/28472523/). To exclude possible drastic structural changes, the pipeline checks the fold change between the synthetic sRNAs and the specified wild-type sRNA using [RNApdist](https://pubmed.ncbi.nlm.nih.gov/22115189/). Both binding energy and distance to wild-type are combined into an comprehensive internal score, on which potential synthetic sRNAs are ranked. To ensure that predicted sRNAs are compatible with possible downstream applications, this ranking is first checked and filtered for illegal sequences such as TypeIIs recognition sites (`exclude_sequence_path` in the config.yml). Subsequently, the ranking is checked for possible off-target effects. For this purpose, SEEDling compares the synthetic sRNAs with the whole target genome via [BLASTn](https://pubmed.ncbi.nlm.nih.gov/20003500/), where the user can determine the threshold via the E-Value parameter. 


## Requirements
The tool was tested on Ubuntu 22.04.1 LTS with the following installations:   
- [Python 3.9.13](https://www.python.org/downloads/release/python-3913/)   
- [conda 4.13.0](https://github.com/conda/conda/releases/tag/4.13.0)    
- [IntaRNA 3.3.1 ](https://github.com/BackofenLab/IntaRNA/releases/tag/v3.3.1)    
- [ViennaRNA 2.4.17](https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.4.17)  
- [NCBI-Blast 2.12.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/)
- [PyYAML 6.0](https://pyyaml.org/)
- [Biopython 1.79](https://biopython.org/wiki/Download)
- [pandas 1.5.1](https://pandas.pydata.org/)


## Converting GenBank
SEEDling relies on GenBank files to function properly. If you have a fasta sequence and a gff annotation file, you can still use SEEDling by converting these files into GenBank format.
One way to accomplish this is by utilizing [EMBOSS seqret](https://www.bioinformatics.nl/cgi-bin/emboss/help/seqret) using the following command.

```
seqret -sequence example_fasta.fasta -feature -fformat gff -fopenfile example_gff.gff -osformat genbank -auto
```
Alternatively, sequence editors like [Benchling](https://www.benchling.com/) and [Geneious Prime](https://www.geneious.com/features/) can be used.

## Usage
## Docker
Due to the number of dependencies, it is recommended to run the pipeline via the Docker Image we've created. 
The SEEDling Docker image can be downloaded from [OwnCloud](https://owncloud.gwdg.de/index.php/s/7tKXsNXfq9OdQzs) or [DockerHub](https://hub.docker.com/r/cedkb/digger_bac-seedling). 
### Running SEEDling
This image was build and tested on Docker version 20.10.23, build 7155243 and Docker Desktop 4.17.0.

**OwnCloud**   
Download SEEDling from [OwnCloud](https://owncloud.gwdg.de/index.php/s/7tKXsNXfq9OdQzs) and load it via the following command
```
docker load -i SEEDling.tar
```
To run the image, use the following command:

```
docker run -it \
--mount type=bind,source=homepath,target=/home/DIGGER/SEEDling/input \
seedling:version1.1
```
Make sure to change the `homepath` to the corresponding folder on the host machine.

**DockerHub**   
Alternatively, the image can be pulled from [DockerHub](https://hub.docker.com/r/cedkb/digger_bac-seedling) using the following command.
``` 
docker pull cedkb/digger_bac-seedling
```
The image can then be exectued using the following command
```
docker run -it \
--mount type=bind,source=homepath,target=/home/DIGGER/SEEDling/input \
cedkb/digger_bac-seedling
```

### Using SEEDling

Change the files in the `input` folder as needed and set SEEDlings parameters in the corresponding config.yml. Note: Make sure that all *paths* in the config.yml include "input/" as a prefix. 

Execute SEEDling by running   
> python3 SEEDling.py -c input/config.yml

SEEDling now executes and generates an output file as defined in the *config.yml*



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
| srna_template | string | GCCACTGCTTTTCTTTGATGTCCCCATTTTG | Template scaffold used for RNApdist |
| exclude_sequences_path | string | input/exclude_sites.fasta | Path to FASTA file that contains sequences to exclude (i.e. restriction recognition sites) |
| include_genes_path | string | input/filter_short.txt | Path to *.txt file of gene names that should be chosen (If left empty all genes will be chosen)| 
| blast_evalue | float | 0.04 | E-value for BLASTn offtarget checks (on target genome) |

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

## License
This is a repository written under the [CC BY-NC-SA 4.0 license](https://creativecommons.org/licenses/by-nc-sa/4.0/)


## References
Martin Mann, Patrick R. Wright, Rolf Backofen, IntaRNA 2.0: enhanced and customizable prediction of RNA–RNA interactions, Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W435–W439, https://doi.org/10.1093/nar/gkx279   

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. In BMC Bioinformatics (Vol. 10, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1471-2105-10-421

Peter J. A. Cock, Tiago Antao, Jeffrey T. Chang, Brad A. Chapman, Cymon J. Cox, Andrew Dalke, Iddo Friedberg, Thomas Hamelryck, Frank Kauff, Bartek Wilczynski, Michiel J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics, Bioinformatics, Volume 25, Issue 11, 1 June 2009, Pages 1422–1423, https://doi.org/10.1093/bioinformatics/btp163

Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. In Algorithms for Molecular Biology (Vol. 6, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1748-7188-6-26
