# Index_picker
This repo is dedicated to a small untility program that help to choose a color-balanced set of index. 

## Background and motivation
Illumina has created a handful of chemistry to generate signal during sequencing. To obtain a high quality base call, the sequencer has to be clearly distinguish cluster on the tile, which can be easily done if the flowcell has high complexity. This assumption may be violated when sequencing low complexity part of the library, such as index. It is therefore important to choose suitable indexes combination that balance in bases such that the complexity is maximized. 

## Idea for implementation
At the moment of writing, there are 2 approaches to solve the problem. The first idea is to brute force the complexity maximization step, and make the best choice by considering the whole list of indexes.
The second idea is to base on the existing list of chosen index, calculate the current base composition, then calculate the potential list of index such that it compensate the base composition and maximize the base balance. 


### Test Case
1. User already have a list of index (dual index, 8bp)
2. User have X more indexes to choose from a defined list

### Usage
Please keep the header of __CD_index.txt__ and __chosen_indexes.txt__.
```
python index_picker.py CD_index.txt chosen_indexes.txt
```

### TODOs
1. provide before and after base composition of index. 
2. print warnings if any of the base exceed limit
