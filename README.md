![logo](doc/images/MetaCHIP_logo.jpg)

Publication
---
+ Song W, Wemheuer B, Zhang S, Steensen K, Thomas T. (2018) MetaCHIP: community-level horizontal gene transfer identification through the combination of best-match and explicit phylogenetic tree approaches (submitted)

+ Contact: Weizhi Song (songwz03@gmail.com), Torsten Thomas(t.thomas@unsw.edu.au)

+ Affiliation: The Centre for Marine Bio-Innovation (CMB), The University of New South Wales, Sydney, Australia

Dependencies:
---

1. For Get_clusters.py
    + Python package:
    [BioPython](https://github.com/biopython/biopython.github.io/),
    [ETE3](http://etetoolkit.org),
    [Numpy](http://www.numpy.org),
    [SciPy](https://www.scipy.org),
    [Matplotlib](http://matplotlib.org)
    + [R](https://www.r-project.org) and its module:
    [optparse](https://cran.r-project.org/web/packages/optparse/index.html),
    [ape](https://cran.r-project.org/web/packages/ape/index.html)
    + [Prodigal](https://github.com/hyattpd/Prodigal)
    + [HMMER](http://hmmer.org)
    + [MAFFT](https://mafft.cbrc.jp/alignment/software/)
    + [FastTree 2.1.9](http://www.microbesonline.org/fasttree/)

1. For Best-match.py
    + Python package:
    [ReportLab](http://www.reportlab.com)
    + [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

1. For Phylogenetic.py
    + Python package:
    [ETE3](http://etetoolkit.org),
    [Pillow](https://pypi.python.org/pypi/Pillow/3.3.1)
    + [Get_Homologues](https://github.com/eead-csic-compbio/get_homologues)
    + [Ranger-DTL 1.0](http://compbio.mit.edu/ranger-dtl/)


How to install:
---
1. Download MetaCHIP scripts and decompress.
1. Isntall required softwares and Python modules.
1. Modify the config file according to the location of required softwares.
1. MetaCHIP is ready to run now.

Notes:
---
1. To get a reliable clusrering results of the input genome bins, their completeness need to be higher than 40% (Figure 3 in doc/MetaCHIP_manuscript.pdf).
1. You can use [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) to get the taxonomy of your input genome bin to further refine the cluster results.
1. The only input for MetaCHIP is a folder which holds the sequence file of a set of genome bins.
1. [Get_Homologues](https://github.com/eead-csic-compbio/get_homologues) is needed to get the Ortholog groups within your input genomes.

How to run:
---

1. Customize file config.txt, specify the path to needed executables in your system.

1. Cluster your input genome bins
        
        python Get_clusters.py -i human_gut_bins -x fasta -p human_gut

1.  You may want to manually modify the grouping results ([prefix]_grouping_g[num].txt) based on the taxonomy classification of your input genome bins

1.  Best-match approach

        python Best-match.py -p human_gut -num_threads 9

1.  Phylogenetic approach

        python Phylogenetic.py -p human_gut


Output files:
---

1. All identified candidates (this file will be used as input for the second step).

        Gene_1	Gene_2	Genome_1_ID	Genome_2_ID	Identity    End_break	Direction
        AAM_03063	BAD_01456	A_1	B_1	100.0	no	AAM<-BAD
        AMAU_05759	BHS_02873	A_4	B_4	79.02	no	AMAU<-BHS
        BAD_00475	AAM_01658	B_1	A_1	74.354	no	BAD<-AAM
        BDS_01877	AMAC_00215	B_2	A_3	100.0	no	BDS<-AMAC
        BGC_01656	AKV_01272	B_3	A_2	100.0	no	BGC<-AKV
        BNM_00983	AMAC_00159	B_5	A_3	75.269	no	BNM<-AMAC
        BNM_02093	AMS_03378	B_5	A_5	100.0	no	BNM<-AMS
        BNM_02445	AMS_01681	B_5	A_5	77.961	no	BNM<-AMS
        BNM_02717	AAM_02737	B_5	A_1	74.47	no	BNM<-AAM

1. Identity distribution plot for each group pair.

    ![identity_distribution](doc/images/identity_distribution.png)

1. ACT image for checking flanking regions of identified HGTs.

    ![flanking_regions](doc/images/flanking_regions.jpg)

1. Combined species tree and gene tree, as well as Ranger-DTL predicted HGTs.

    ![Combined_tree](doc/images/Combined_trees.png)
