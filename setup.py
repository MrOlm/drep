from setuptools import setup, find_packages
import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'drep', 'VERSION'))
    return versionFile.read().strip()

setup(name='drep',
      version=version(),
      description='De-replication of microbial genomes assembled from multiple samples',
      url='https://github.com/MrOlm/drep',
      author='Matt Olm',
      author_email='mattolm@berkeley.edu',
      license='MIT',
      package_data={'drep': ['VERSION']},
      packages=find_packages(),
      scripts=['bin/dRep', 'helper_scripts/parse_stb.py', 'helper_scripts/ScaffoldLevel_dRep.py'],
      install_requires=[
          'numpy',
          'pandas',
          'seaborn',
          'matplotlib',
          'biopython',
          'scikit-learn',
          'tqdm',
          'pytest'
      ],
      zip_safe=False)
