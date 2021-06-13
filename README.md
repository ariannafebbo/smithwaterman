# SmithWaterman
A pairwise local sequence alignment method, to find the optimal local sequence alignment between two nucleotide sequences, using the Smith-Waterman algorithm. 
The **SmithWaterman** project is available at [https://github.com/ariannafebbo/smithwaterman]

## Description
The outcome of the **SmithWaterman** project is a function which takes as parameters: 
1. The first nucleotide sequence `seq1`;
2. the second nucleotide sequence `seq2`;
3. the cost of a match `match_cost`;
4. the cost of a mismatch `mismatch_cost`;
5. the cost of a gap/indel `gap_cost`;

and returns as result the optimal alignment of the two input sequences and the alignment score.

## How to install the package
Once having downloaded the package **SmithWaterman**, type from terminal, inside the folder **SmithWaterman**:
- `pip install .`

Then import the package in your python script and type:
- `from smith_watermann import sw`

### Usage example
- `sw.find_alignment("ACTG","AGTG")`

`match_cost`, `mismatch_cost` and `gap_cost` are respectevely set as 1, -1 , -2 as default parameters.

## How to run test
Type from terminal, inside the folder **SmithWaterman**:
- `pip install pipenv` builds the ad hoc environment;
- `pipenv shell` install dependencies;
- `pipenv update` updates the environment;
- `pytest` run tests.






