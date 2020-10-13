from setuptools import setup, find_packages

setup(name='sniphles',
      version='0.1',
      description='''sniphles is a simple python package for .''',
      url='https://github.com/collaborativebioinformatics/Sniphles',
      author='The Pan-Structural Variation hackathon',
      author_email='decosterwouter@gmail.com',
      license='MIT',
      ### Dependencies
      install_requires=[
          'pysam',
          'numpy',
      ],
      packages=find_packages(),
      python_requires='>=3.8',
      zip_safe=False)
