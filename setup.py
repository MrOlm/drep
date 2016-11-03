from setuptools import setup

setup(name='drep',
      version='0.1.1',
      description='De-replication of microbial genomes assembled from multiple time-points',
      url='https://github.com/MrOlm/drep',
      author='Matt Olm',
      author_email='mattolm@berkeley.edu',
      license='MIT',
      packages=['drep'],
      scripts=['bin/Drep'],
      install_requires=[
          #'numpy',
          #'pandas',
      ],
      zip_safe=False)
