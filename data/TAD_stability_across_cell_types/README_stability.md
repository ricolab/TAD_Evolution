The stability of a TAD corresponds to the number of cell types where its boundaries are present. To study the stability of TADs across cell types, we used the TAD maps and associated stability scores calculated by McArthur and Capra (https://pubmed.ncbi.nlm.nih.gov/33545030/ and https://zenodo.org/records/4156731). Their dataset includes definitions for 14,345 TAD boundaries and assigns to each of them a parameter, the stability percentile, that is associated to the number of cell-types where it is present (which varies between one and 37). For example, if a TAD boundary is present in only one cell-type, it corresponds to a stability percentile of 0.113175, but if it is present in 24 cell-types, then its stability percentile is 0.901778. All TAD boundaries in this dataset span for 100 kb.

We separated the dataset of TAD boundaries in two groups of approximately the same size depending on their presence in different cell-types:

- the "stable" group, which was made of TAD boundaries present in five or more different cell-types, including 7,509 TAD boundaries (52%).
- and the "unstable" group, which was made using the remaining 6,836 TAD boundaries (48%), present in less than five different cell-types.

For each of these two groups, we created a list of TADs as the region between two consecutive TAD boundaries, by taking the central value between the beginning and end of each TAD boundary definition. In addition, we removed all the TADs that were nesting a TAD boundary from the other group (so that the final TAD definitions in either group were empty of other TAD boundaries). Using these criteria, we obtained two datasets:

- one with 4749 stable TADs, with an average length of 154 kb.
- another one with 4086 unstable TADs, with an average length of 268 kb.
