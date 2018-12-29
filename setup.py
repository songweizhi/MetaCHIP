from setuptools import setup, find_packages


__long_description__ = '''

MetaCHIP: community-level HGT identification pipeline

Weizhi Song (songwz03@gmail.com)

The Centre for Marine Bio-Innovation (CMB), 
University of New South Wales, Sydney, Australia

'''

__version__ = '1.0.1'

setup(name="MetaCHIP",
      version=__version__,
      long_description=__long_description__,
      license="GPL3+",
      author="Weizhi Song, Torsten Thomas",
      author_email="songwz03@gmail.com",
      keywords="Bioinformatics Metagenomics HGT_detection",
      description="HGT detection pipeline",
      url="https://github.com/songweizhi/MetaCHIP",
      packages=find_packages(exclude=['contrib','docs']),
      package_data={'': ['*.r', '*.R', '*.py', '*.hmm']},
      include_package_data=True,
      install_requires=['biopython', 'matplotlib', 'numpy', 'scipy', 'reportlab', 'ete3'],

      scripts=['bin/MetaCHIP',
               'MetaCHIP/PI.py',
               'MetaCHIP/BM.py',
               'MetaCHIP/PG.py',
               'MetaCHIP/MetaCHIP_config.py',
               'MetaCHIP/MetaCHIP_phylo.hmm',
               'MetaCHIP/MetaCHIP_circos_HGT.R',
               'MetaCHIP/MetaCHIP_add_group_to_tree.R'])

