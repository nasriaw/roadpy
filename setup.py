#setup.py
from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'roadpy'
LONG_DESCRIPTION = 'package for roadway geometric design, '

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="roadpy", 
        version=VERSION,
        author="nasri",
        author_email="<nasriaw@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=(['mathematical'], ['pandas'], ['haversine'], ['numpy'], ['matplotlib']), # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'roadpy'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
      )