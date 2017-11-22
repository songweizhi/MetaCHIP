# Steps to install ete3 

# 1. Install PyQt and Sip through Anaconda (https://repo.continuum.io/archive/Anaconda3-4.2.0-MacOSX-x86_64.sh)

# 2. Install Anaconda
bash Anaconda3-4.2.0-MacOSX-x86_64.sh

# 3. add conda to path
export PATH=~/anaconda/bin:$PATH

# 4. install the right version of pyqt (https://github.com/etetoolkit/ete/issues/195)
conda install pyqt=4

# 5. remove python2 from PATH:
export PATH="/Users/songweizhi/anaconda/bin:/Users/songweizhi/anaconda3/bin:/Users/songweizhi/miniconda3/bin:/Library/Frameworks/Python.framework/Versions/3.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/ncbi/blast/bin"

# 6. install ete3 through conda
conda install -c etetoolkit ete3 ete3_external_apps

# 7. or install ete3 through pip3
pip3 install --upgrade ete3




# Customize NCBI taxonomy database

# 1. Get customized taxdump.tar.gz
tar -czf taxdump.tar.gz names.dmp nodes.dmp merged.dmp 

# 2. modify ncbiquery.py (/Users/songweizhi/anaconda3/lib/python3.5/site-packages/ete3/ncbi_taxonomy/ncbiquery.py)
Delete line 720: urlretrieve("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", "taxdump.tar.gz")
Change the path in line 722 to your customized "taxdump.tar.gz"
