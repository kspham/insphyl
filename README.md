#insphyl
##Introduction
Construct insertion-based phylogenetic tree
##Tool Description
Here are some parameters for `insphyl`  
* `-ref`	\[reference genome\]	The reference file in FASTA format. The file's name must have the tail of `.fasta`. This parameter is mandatory.
* `-assembly` \[assembly files\]	The assembly file(s) of the sample. Each sample must be put in a single FASTA file. The file's name must have the tail of `.fasta`. This parameter is mandatory.
* `-output`	(output directory)	This is where all output files of `insphyl` will be created. Default is `./`.
* `-minlen`	(minimum length)	The minimum length of the insertions. Default is `500`.
* `-id` (minimum identity)	The minimum identity for the clustering process. Default is `0.9`.  

Example
```
insphyl -ref reference.fasta -assembly sample_1.fasta sample_2.fasta sample_3.fasta -minlen 700 -id 0.95 -output project_1/analysis
```
##Pipeline
1. Get ID from ENA - `qualify_assembly.py`  
2. Download rawdata - `download_rawdata.py`  
3. Trim - `trim.py`  
4. Assemble - `assembly.py`   
5. Qualify Assembly - `qualify_assembly.py`  
6. Get NSR - `get_NSR.py`  
7. Format NSR - `format_NSR.py`  
8. Cluster NSR - `cluster_NSR.py`  
9. Cluster for all isolates - `cluster_all.py`  
10. Index clusters - `index_cluster.py`  
11. Calculate distance in array - `distance_array.py`  
12. Convert to MEG for neighbor-joining - `convert_to_MEG.py`  
13. Remap with rawdata - `remap.py`  
14. Construct heatmap of converage for all isolates - `coverage_heatmap.py`
15. Pick out epidemic insertions - `epidemic.py`
16. Annotate True Possitive result - `annotate.py`

##Process
__Data Preparation__ 1-4  
__Tool__ 6-12  
__Analysis__ 11-16  
##Script Description
###`get_ID_from_ENA.py`
Collect IDs and ftp links. Output file is a tab-delimited file of `id_list.txt`
###`download_rawdata.py`
Download all sequencing raw data based on the links provided in the `id_list.txt`
###`trim.py`
Trim the raw data with `trimmomatic` with the following setting:
```
PE  
-phred33  
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10  
LEADING:3  
TRAILING:3  
SLIDINGWINDOW:4:15  
MINLEN:36
```
###`assembly.py`
A script for genome assembly of pair-end SRA based on `SPAdes`
###`qualify_assembly.py`
Generate statistical report on assemblies by using `QUAST`.
###`get_NSR.py`
Index the systenic region pairwisely for all isolate with the reference of `HE681097.fasta`, an EMRSA-15. Then reverse all index to get NSRs as FASTA and at the same time, generate all variants in VCF using mapping of assemblies.
###`format_NSR.py`
Rename all NSG respective to the experimental name with id starting from 1. Part of this script is also used in clustering scripts.
###`cluster_NSR.py`
To get rid of the multi-resolution of C-Sibelia output, this script is used for remove all NSRs that support for the same inserted DNA segment. The threshold of identity is set at 90
###`cluster_all.py`
Cluster all insertions from all isolate into a single FASTA of `all.fasta`. The clusters in this file are not merged into centroid sequences. This output is imperative for distance calculation.
###`index_cluster.py`
Format all heads of the cluster into a combination of their sequences. The Id will also be set from 1.
###`distance_array.py`
From indexed FASTA file of the clusters, this script will generate an array represent for all isolates and their insertions. The value for absence and presence of a variant is set at 0 and 1 respectively.
###`convert_to_MEG.py`
This script is for the conversion of distance array to MEG format of lower left array.
###`remap.py`
Map all trimmed SRA to the insertions. The script uses BWA and Samtools to generate BAM file.
###`coverage_heatmap.py`
Create an array that records the breadth of coverage from the remapped BAM file. The script is assisted by `bedtools`.
###`epidemic.py`
Based on the record of mapping coverage, the script pick out the insertions that shared only by the outbreak isolates. The sequences is saved in FASTA format.
###`annotate.py`
The script uses `prokka` to annotate the epidemic insertions
##Author
**Tri Le**  
Intenational University, Ho Chi Minh, Vietnam  
Institute of Computer Science and Technology, Ho Chi Minh, Vietnam  
tri.lqm@icst.org.vn
