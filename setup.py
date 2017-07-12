from setuptools import setup
import Pandoo
import os

def read(fname):
    '''
    Read the README
    '''
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'pandoo',
    version = Pandoo.__version__,
    description = Pandoo.__description__,
    long_description=read('README'),
    classifiers = ['Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: GNU Affero General ' +
                   'Public License v3 or later (AGPLv3+)',
                   'Programming Language :: Python :: 3.5',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Medical Science Apps.',
                   'Intended Audience :: Science/Research'],
    keywords = ['pipeline',
                'ruffus',
                'bacteria',
                'contigs',
                'assembly',
                'assemblies',
                'reads',
                'paired-end',
                'public health microbiology',
                'microbial genomics'],
    download_url = Pandoo.__download_url__,
    author = Pandoo.__author__,
    author_email = Pandoo.__author_email__,
    license = Pandoo.__license__,
    packages = ['Pandoo'],
    scripts = ['Pandoo/pandoo'],
    include_package_data = True,
    install_requires = ['scipy==0.19.1',
                        'numpy==1.13.1',
                        'pandas==0.20.3',
                        'ete3==3.0.0b35',
                        'ruffus==2.6.3']
    )