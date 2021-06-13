from setuptools import setup

setup(
    # Needed to silence warnings
    name='SmithWaterman',
    url='https://github.com/ariannafebbo/smithwaterman',
    author='Arianna Febbo',
    # Needed to actually package something
    packages=['smith_waterman'],
    # Needed for dependencies
    install_requires=['numpy'],
    # *strongly* suggested for sharing
    version='0.5',
    license='MIT',
    description='A pairwise local sequence alignment method, to find the optimal global sequence alignment between two nucleotide sequences, using the Smith-Waterman algorithm',
    long_description=open('README.md').read(),
)