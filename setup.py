from setuptools import setup

setup(name='drep',
      version='0.2.0',
      description='De-replication of microbial genomes assembled from multiple samples',
      url='https://github.com/MrOlm/drep',
      author='Matt Olm',
      author_email='mattolm@berkeley.edu',
      license='MIT',
      packages=['drep'],
      scripts=['bin/dRep'],
      install_requires=[
          #'numpy',
          #'pandas',
      ],
      zip_safe=False)
