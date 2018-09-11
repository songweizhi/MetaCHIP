![logo](images/MetaCHIP_logo.jpg)


Publication
---
+ Song WZ, Wemheuer B, Zhang S, Steensen K, Thomas T* (2018) MetaCHIP: community-level horizontal gene transfer identification through the combination of best-match and explicit phylogenetic tree approaches. Microbiome **(Under review)**.
+ Contact: Weizhi Song (songwz03@gmail.com), Torsten Thomas(t.thomas@unsw.edu.au)

+ Affiliation: The Centre for Marine Bio-Innovation (CMB), The University of New South Wales, Sydney, Australia


Dependencies:
---

+ Python packages:
[BioPython](https://github.com/biopython/biopython.github.io/),
[ETE3](http://etetoolkit.org),
[Numpy](http://www.numpy.org),
[SciPy](https://www.scipy.org),
[Matplotlib](http://matplotlib.org)
[ReportLab](http://www.reportlab.com)

+ R packages:
[optparse](https://cran.r-project.org/web/packages/optparse/index.html),
[ape](https://cran.r-project.org/web/packages/ape/index.html)
[circlize](https://cran.r-project.org/web/packages/circlize/index.html)

+ [Prodigal](https://github.com/hyattpd/Prodigal)
+ [HMMER](http://hmmer.org)
+ [MAFFT](https://mafft.cbrc.jp/alignment/software/)
+ [FastTree](http://www.microbesonline.org/fasttree/)
+ [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
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
1. To get a reliable clusrering results of the input genome bins, their completeness need to be higher than 40% (Figure 3 in manual/MetaCHIP_manuscript.pdf).
1. You can use [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) to get the taxonomy of your input genome bin to further refine the cluster results.
1. The only input for MetaCHIP is a folder which holds the sequence file of a set of genome bins.
1. [Get_Homologues](https://github.com/eead-csic-compbio/get_homologues) is needed to get the Ortholog groups within your input genomes.


How to run:
---
+ A detailed manual can be found at [manual/MetaCHIP_User_Manual.pdf](https://github.com/songweizhi/MetaCHIP/blob/master/manual/MetaCHIP_User_Manual.pdf).


Output files:
---
1. Identified HGT candidates.

    |Gene_1|Gene_2|Genome_1_ID|Genome_2_ID|Identity|End_break|Direction|
    |---|---|---|---|---|---|---|
    |AAM_03063|BAD_01456|A_1|B_1|100.0|no|AAM<-BAD|
    |AMAU_05759|BHS_02873|A_4|B_4|79.02|no|AMAU<-BHS|
    |BAD_00475|AAM_01658|B_1|A_1|74.354|no|BAD<-AAM|
    |BDS_01877|AMAC_00215|B_2|A_3|100.0|no|BDS<-AMAC|
    |BGC_01656|AKV_01272|B_3|A_2|100.0|no|BGC<-AKV|
    |BNM_00983|AMAC_00159|B_5|A_3|75.269|no|BNM<-AMAC|
    |BNM_02093|AMS_03378|B_5|A_5|100.0|no|BNM<-AMS|
    |BNM_02445|AMS_01681|B_5|A_5|77.961|no|BNM<-AMS|
    |BNM_02717|AAM_02737|B_5|A_1|74.47|no|BNM<-AAM|

1. Identity distribution for each group pair.

    ![identity_distribution](images/identity_distribution.png)

1. Flanking regions of identified HGTs.

    ![flanking_regions](images/flanking_regions.jpg)

1. Gene flow between groups.

    ![Gene_flow](images/Gene_flow.jpg)

