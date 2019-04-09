import os
from setuptools import setup, find_packages


def version():

    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'MetaCHIP', 'VERSION'))

    return version_file.readline().strip()


__long_description__ = '''

MetaCHIP: community-level HGT identification pipeline

Weizhi Song (songwz03@gmail.com)

The Centre for Marine Bio-Innovation (CMB), 
University of New South Wales, Sydney, Australia

'''


setup(name="MetaCHIP",
      version=version(),
      long_description=__long_description__,
      license="GPL3+",
      author="Weizhi Song, Torsten Thomas",
      author_email="songwz03@gmail.com",
      keywords="Bioinformatics Metagenomics HGT_detection",
      description="HGT detection pipeline",
      url="https://github.com/songweizhi/MetaCHIP",
      packages=['MetaCHIP'],
      package_data={'': ['*.r', '*.R', '*.py', '*.hmm', 'VERSION', 'Ranger-DTL-Dated.mac', 'Ranger-DTL-Dated.linux']},
      include_package_data=True,
      install_requires=['biopython', 'matplotlib', 'numpy', 'scipy', 'reportlab', 'ete3'],
      scripts=['bin/MetaCHIP'])

