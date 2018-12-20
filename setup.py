from setuptools import setup, find_packages


setup(
    name="MetaCHIP",
    version="1.0",
    license="GPL3+",
    author="Weizhi Song, Torsten Thomas",
    author_email="songwz03@gmail.com",
    keywords="Bioinformatics Metagenomics HGT",
    description="HGT detection pipeline",
    url="https://github.com/songweizhi/MetaCHIP",



    packages=find_packages(exclude=['contrib','docs']),

    scripts=['MetaCHIP.py',
             'scripts/PI.py',
             'scripts/BM.py',
             'scripts/PG.py'],

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['biopython',
                      'matplotlib',
                      'numpy',
                      'scipy',
                      'reportlab',
                      'ete3'],

    package_data={
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.txt', '*.rst'],
        # And include any *.msg files found in the 'hello' package, too:
        'hello': ['*.msg'],
    },

    # metadata to display on PyPI

    project_urls={
        "Bug Tracker": "https://bugs.example.com/HelloWorld/",
        "Documentation": "https://docs.example.com/HelloWorld/",
        "Source Code": "https://code.example.com/HelloWorld/",
    }

    # could also include long_description, download_url, classifiers, etc.
)
