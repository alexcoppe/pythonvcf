import setuptools


# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setuptools.setup(
     name='pythonvcf.py',  
     version='0.0.1',
     author="Alessandro Coppe",
     author_email="",
     description="A python package for manipulatin VCF, variant call format, files",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="",
     scripts=["scripts/pythonvcf.py"],
     install_requires=[
      ],
     classifiers=[
         "Programming Language :: Python :: 3",
         "Operating System :: OS Independent",
     ],
)
