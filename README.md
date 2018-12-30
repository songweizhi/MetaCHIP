![logo](images/MetaCHIP_logo.jpg)


Publication
---

+ **Song WZ**, Wemheuer B, Zhang S, Steensen K, Thomas T* (2018) MetaCHIP: community-level horizontal gene transfer identification through the combination of best-match and explicit phylogenetic tree approaches. Microbiome **(Under review)** [PDF](https://songweizhi.github.io/assets/pdfs/Publication_2018_MetaCHIP_manuscript.pdf) 
+ Contact: Weizhi Song (songwz03@gmail.com), Torsten Thomas (t.thomas@unsw.edu.au)
+ Affiliation: The Centre for Marine Bio-Innovation, University of New South Wales, Sydney, Australia


Dependencies:
---

### Python libraries
* [BioPython](https://github.com/biopython/biopython.github.io/): Python tools for computational molecular biology.
* [Numpy](http://www.numpy.org): fundamental package for scientific computing with Python.
* [SciPy](https://www.scipy.org): Python-based ecosystem for mathematics, science, and engineering.
* [Matplotlib](http://matplotlib.org): Python plotting library.
* [ETE3](http://etetoolkit.org): Python environment for tree exploration.

### R packages
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html): command line option parser in R.
* [ape](https://cran.r-project.org/web/packages/ape/index.html): package for analyses of phylogenetics and evolution in R.
* [circlize](https://cran.r-project.org/web/packages/circlize/index.html): package for circular visualization.

### Third-party software
MetaCHIP makes use of the following 3rd party dependencies and assumes these are on your system path. Specify full path 
to their executables in the config file if they are not on the system path.  
* [Prodigal](https://github.com/hyattpd/Prodigal): protein-coding gene prediction tool for prokaryotic genomes.
* [HMMER](http://hmmer.org): tool for biosequence analysis using profile hidden Markov models.
* [MAFFT](https://mafft.cbrc.jp/alignment/software/): multiple sequences alignment program.
* [FastTree](http://www.microbesonline.org/fasttree/): tool for inferring phylogenies from alignments .
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download): you know what it is!
* [Ranger-DTL 2.0](https://compbio.engr.uconn.edu/software/RANGER-DTL/): software for inferring gene family evolution.


How to install:
---

1. MetaCHIP can be installed with pip after you have its 3rd party dependencies and R packages installed on your system. 
Python libraries required by MetaCHIP will be installed automatically during the pip installation:

        pip install MetaCHIP
        
1. You can either add MetaCHIP's dependencies to your system path or specify full path to their executables in MetaCHIP_config.py, which is in folder lib/site-packages/MetaCHIP.


How to run:
---

1. A detailed manual can be found at [manual/MetaCHIP_User_Manual_v1.0.1.pdf](https://github.com/songweizhi/MetaCHIP/blob/master/manual/MetaCHIP_User_Manual_v1.0.1.pdf).

1. The input files for MetaCHIP include a folder that holds the sequence file [[example](https://github.com/songweizhi/MetaCHIP/blob/master/input_file_examples/human_gut_bins)] 
of all query genomes, as well as a text file [[example](https://github.com/songweizhi/MetaCHIP/blob/master/input_file_examples/human_gut_bins_GTDB.tsv)], 
which holds taxonomic classification of all input genomes. Please make sure the length of sequence IDs for all input genomes is **NO LONGER THAN** 22 letters.

1. [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) is recommended for taxonomic classification of input genomes.

1. Options for argument '-r' in module PI, BM and PG include: d (domain), p (phylum), c (class), o (order), f (family) and g (genus).

1. Output format for BLASTN in the PI and BM steps: 
        
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

1. Some examples: 

        # show help information
        $ MetaCHIP -h
        $ MetaCHIP PI -h
        $ MetaCHIP BM -h
        $ MetaCHIP PG -h
        
        # run MetaCHIP by grouping input genomes at Class level according to their taxonomic classifications
        $ MetaCHIP PI -i bin_folder -x fasta -taxon NorthSea_GTDB_output.tsv -r c -p NorthSea -t 6
        $ MetaCHIP BM -p NorthSea -r c -t 6
        $ MetaCHIP PG -p NorthSea -r c -t 6

        # run MetaCHIP with customized grouping profile
        $ MetaCHIP PI -i bin_folder -x fasta -g customized_grouping.txt -p NorthSea -t 6
        $ MetaCHIP BM -p NorthSea -g customized_grouping.txt -t 6
        $ MetaCHIP PG -p NorthSea -g customized_grouping.txt -t 6
        
        # Detect HGT with MetaCHIP at multiple levels (e.g. class and order)
        $ MetaCHIP PI -i bin_folder -x fasta -taxon NorthSea_GTDB_output.tsv -r c -p NorthSea -t 6
        $ MetaCHIP PI -i bin_folder -x fasta -taxon NorthSea_GTDB_output.tsv -r o -p NorthSea -t 6 -grouping_only
        $ MetaCHIP BM -p NorthSea -r c -t 6
        $ MetaCHIP BM -p NorthSea -r o -t 6
        $ MetaCHIP PG -p NorthSea -r c -t 6
        $ MetaCHIP PG -p NorthSea -r o -t 6


Output files:
---

1. Identified HGT candidates and their sequences. See some examples of contig end_match and full_length_match below. 

    |Gene_1|Gene_2|Gene_1_group|Gene_2_group|Identity|end_match|full_length_match|Direction|
    |---|---|---|---|---|---|---|---|
    |AAM_03063|BAD_01456|A|B|100.0|no|no|AAM<-BAD|
    |AMAU_05759|BHS_02873|A|B|79.02|yes|no|AMAU<-BHS|
    |BAD_00475|AAM_01658|B|A|74.354|no|no|BAD<-AAM|
    |BDS_01877|AMAC_00215|B|A|100.0|no|yes|BDS<-AMAC|
    |BGC_01656|AKV_01272|B|A|100.0|no|no|BGC<-AKV|
    |BNM_00983|AMAC_00159|B|A|75.269|no|no|BNM<-AMAC|
    |BNM_02093|AMS_03378|B|A|100.0|yes|no|BNM<-AMS|
    |BNM_02445|AMS_01681|B|A|77.961|no|yes|BNM<-AMS|
    |BNM_02717|AAM_02737|B|A|74.47|no|no|BNM<-AAM|

1. Flanking regions of identified HGT candidates.
    ![flanking_regions](images/flanking_regions.jpg)
        
1. Gene flow between groups.
    ![Gene_flow](images/Gene_flow.jpg)
    
1. Identity distribution of identified HGTs.
    ![HGT_identity_distribution](images/HGT_identity_distribution.png)

1. The number of identified HGT in each genome/group.
    ![HGT_per_genome_group](images/HGT_per_genome_group.png)
 
1. Examples of contig end match.
    ![end_match_1](images/end_match_1.jpg)   
    ![end_match_2](images/end_match_2.jpg)
        
1. Examples of full-length contig match
    ![full_length_match_1](images/full_length_match_1.jpg)
    ![full_length_match_2](images/full_length_match_2.jpg)
    