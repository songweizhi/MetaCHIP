![logo](images/MetaCHIP_logo.jpg)


Publication:
---

+ **Song WZ**, Wemheuer B, Zhang S, Steensen K, Thomas T* (2019) MetaCHIP: community-level horizontal gene transfer identification through the combination of best-match and phylogenetic approaches. Microbiome. 7:36 https://doi.org/10.1186/s40168-019-0649-y
+ Contact: Weizhi Song (songwz03@gmail.com), Torsten Thomas (t.thomas@unsw.edu.au)
+ Affiliation: The Centre for Marine Bio-Innovation, University of New South Wales, Sydney, Australia


Dependencies:
---

### Python libraries (will be installed automatically during pip installation)
* [BioPython](https://github.com/biopython/biopython.github.io/): Python tools for computational molecular biology.
* [Numpy](http://www.numpy.org): fundamental package for scientific computing with Python.
* [SciPy](https://www.scipy.org): Python-based ecosystem for mathematics, science, and engineering.
* [Matplotlib](http://matplotlib.org): Python plotting library.
* [ETE3](http://etetoolkit.org): Python environment for tree exploration.

### R packages (will be installed automatically when needed)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html): command line option parser in R.
* [ape](https://cran.r-project.org/web/packages/ape/index.html): package for analyses of phylogenetics and evolution in R.
* [circlize](https://cran.r-project.org/web/packages/circlize/index.html): package for circular visualization.

### Third-party software
MetaCHIP makes use of the following 3rd party dependencies and assumes these are on your system path. Specify full path 
to their executables in the config file if they are not on the system path.  
* [Prodigal](https://github.com/hyattpd/Prodigal): protein-coding gene prediction tool for prokaryotic genomes.
* [HMMER](http://hmmer.org): tool for biosequence analysis using profile hidden Markov models.
* [MAFFT](https://mafft.cbrc.jp/alignment/software/): multiple sequences alignment program.
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): you know what it is!
* [FastTree](http://www.microbesonline.org/fasttree/): tool for inferring phylogenies from alignments (part of MetaCHIP, **no need to install**).
* [Ranger-DTL 2.0](https://compbio.engr.uconn.edu/software/RANGER-DTL/): software for inferring gene family evolution (part of MetaCHIP, **no need to install**).


How to install:
---

1. MetaCHIP can be installed with pip after you have its 3rd party dependencies and R packages installed on your system. 
Python libraries required by MetaCHIP will be installed automatically during the pip installation:

        pip install MetaCHIP
        
1. You can either add MetaCHIP's dependencies to your system path or specify full path to their executables in MetaCHIP_config.py which can be found in Python's folder lib/site-packages/MetaCHIP.


How to run:
---

1. The input files for MetaCHIP include a folder that holds the sequence file [[example](https://github.com/songweizhi/MetaCHIP/blob/master/input_file_examples/human_gut_bins)] 
of all query genomes, as well as a text file [[example](https://github.com/songweizhi/MetaCHIP/blob/master/input_file_examples/human_gut_bins_GTDB.tsv)], 
which holds taxonomic classification of all input genomes. Please make sure the length of **sequence ID** for sequences of all input genomes are <font color="red"> **NO LONGER THAN 22 letters** </font>.

1. [**GTDB-Tk**](https://github.com/Ecogenomics/GTDBTk) is recommended for taxonomic classification of input genomes. 
GTDB-Tk with produce two files ([prefix].bac120.summary.tsv and [prefix].ar122.summary.tsv) containing the classification results 
if you have both bacterial and archaeal genomes included in your queries. You need to combine the two files into one and feed it as an input for MetaCHIP.

1. Options for argument '-r' in module PI, BM and PG include: d (domain), p (phylum), c (class), o (order), f (family) and g (genus).

1. Output format for BLASTN in the PI and BM steps: 
        
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

1. Some examples: 

        # show help information
        MetaCHIP -h
        MetaCHIP PI -h
        MetaCHIP BP -h
                
        # Example 1: run MetaCHIP at single level (class)
        MetaCHIP PI -i bin_folder -x fasta -taxon NorthSea_GTDB_bac120_ar122_combined.tsv -r c -p NorthSea -t 6
        MetaCHIP BP -p NorthSea -r c -t 6

        # Example 2: run MetaCHIP at multiple levels (e.g. class and order)
        MetaCHIP PI -i bin_folder -x fasta -taxon NorthSea_GTDB_bac120_ar122_combined.tsv -r co -p NorthSea -t 6
        MetaCHIP BP -p NorthSea -r co -t 6

        # Example 3: run MetaCHIP at multiple levels (e.g. phylum, class, order and family)
        MetaCHIP PI -i bin_folder -x fasta -taxon NorthSea_GTDB_bac120_ar122_combined.tsv -r pcofg -p NorthSea -t 6
        MetaCHIP BP -p NorthSea -r pcofg -t 6

        # Example 4: run MetaCHIP with customized grouping profile
        MetaCHIP PI -i bin_folder -x fasta -g customized_grouping.txt -p NorthSea -t 6
        MetaCHIP BP -p NorthSea -g customized_grouping.txt -t 6
        

Output files:
---

1. A Tab delimited text file containing all identified HGTs. Filename format: [prefix]_[taxon_ranks]_detected_HGTs.txt

    |Column|Description|
    |---|---|
    |Gene_1|The 1st gene involved in a HGT event|
    |Gene_2|The 2nd gene involved in a HGT event|
    |Identity|Identity between Gene_1 and Gene_2|
    |occurence(taxon_ranks)|Only for multiple-level detections. If you performed HGT detection at phylum, class and order levels, a number of "011" means current HGT was identified at class and order levels, but not phylum level.|
    |end_match|End match or not (see examples below)|
    |full_length_match|Full length match or not (see examples below)|
    |direction|Direction of gene flow|   


1. Nucleotide and amino acid sequences of identified donor and recipient genes.


1. Flanking regions of identified HGTs. Genes encoded on the forward strand are displayed in light blue, and genes coded on the reverse strand are displayed in light green. The name of genes predicted to be HGT are highlighted in blue, large font with pairwise identity given in parentheses. Contig names are provided at the left bottom of the sequence tracks and numbers following the contig name refer to the distances between the gene subject to HGT and either the left or right end of the contig. Red bars show similarities of the matched regions between the contigs based on BLASTN results.
    ![flanking_regions](images/flanking_regions.png)

        
1. Gene flow between groups. Bands connect donors and recipients, with the width of the band correlating to the number of HGTs and **the colour corresponding to the donors**.
    ![Gene_flow](images/Gene_flow.jpg)


1. Examples of contig end match.
    ![end_match_1](images/end_match_1.jpg)   
    ![end_match_2](images/end_match_2.jpg)

        
1. Examples of full-length contig match
    ![full_length_match_1](images/full_length_match_1.jpg)
    ![full_length_match_2](images/full_length_match_2.jpg)
