# Index_picker
This repo is dedicated to a small untility program that help choose a color-balanced set of index. 

## Background and motivation
Illumina has created a handful of chemistry to generate signal during sequencing. To obtain a high quality base call, the sequencer has to be clearly distinguish cluster on the tile, which can be easily done if the flowcell has high complexity. This assumption may be violated when sequencing low complexity part of the library, such as index. It is therefore important to choose suitable indexes combination that balance in bases such that the complexity is maximized. 

## Idea for implementation
At the moment of writing, there are 2 approaches to solve the problem. The first idea is to brute force the step maximization step, and make the best choice by considering the whole list of indexes.
The second idea is to base on the existing list of chosen index, calculate the current base composition, then calculate the potential list of index such that it compensate the base composition and maximize the base balance. 

This repo will attempt the second idea, just because it seems more fun to build. While it has not escape my attention that, given the limited number of indexes, the first idea maybe more suitable. It also make sense as it give the global maximum (best) answer to the question. 

### Test Case
1. User already have a list of index (dual index, 8bp)
2. User have X more indexes to choose from a defined list