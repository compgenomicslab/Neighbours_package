# Neighbours_package
**Neighbours_package v2**

This is the new package for neighbours analysis.
The package contains two main scripts:

A) neigh.py  <br />
B) neigh_graphic.py

Both scripts use mongo_client.py to set up mongo connection, and also de folder data, in which is storage KEGG pathways information.

Finally, neigh_functions.py is a module required for neigh.py. Contains most functions used by neigh.py script.

<br />

#### neigh.py usage: ####

```
Usage:  neigh.py -u <unigene cluster> or neigh.py -f <file(unigenes_list).txt> 
```

use `-u` just to inspect one unigene. If you want to analyze a complete list of unigene use `-f` flag

<br />

#### neigh.py OUTPUT ####
The output contains three different lines:

1.  Keegs assignation
2.  EGGnog assignation
3. List of neighbours genes for analysed unigene

<br />

In the first and second lines, columns from 1 to 11 contain:

- **gmgc**: unigene cluster

- **query_cogs**: cogs(Kegg,Eggnog) assigned in mongo database to query unigene cluster 

- **subject_cogs**:  cogs(Kegg,Eggnog) predicted by neighbourhood analysis 

- **analysed_orfs**: number of ORFs that contain the query unigene cluster

- **number_neigh_genes**: the total sum of neighbours genes presents in all the ORFs that conforms the unigene cluster 

- **number_neigh_with_cogs**: number of neigh genes with cog assignation(Kegg or EGGnog) 

- **unique_cogs**: number of unique cogs present in the query unigene cluster 

- **count_of_cogs**: the total sum of all cogs(kegg or EGGnog) presentin the query unigene cluster 

- **cog_conservation**: 1 - [(unique_cogs/count_of_cogs)/number_neigh_with_cogs] 

- **hit_cog_percentage**: percentage of presence of a determined cog in the neighbourhood. 
  cog_description: Functional description of cogs

  <br />
  <br />

For the 3rd line, List of neighbours genes for analysed unigene. The ORFs are separeted by comma, and neighbourhood of every ORF is represented by its unigene identifier and not for its gen ID. For example, Below you can see one selected ORF from the 000_000_005 unigene cluster. Five unigenes are separated by '@', 2 upstream of the ORFs and 2 downstream of the query ORF, and in the middle its present the query unigene that represent the ORF.

`00_334_750@000_472_881@000_000_005@000_472_880@000_334_749`

In the case there are not neigh genes presents, the unigene identifier would be substitute by 'NA'.

`NA@NA@000_000_005@NA@NA`

<br />

### neigh_graphic.py usage: ###

```
Usage: neigh_graphic.py <number_of_neighbours_genes_to_display> <unigene_cluster>
Example: neigh_graphic.py 2 000_000_005
```
<br />

`number_of_neighbours_genes_to_display:`  It controls the number of unigenes to display upstream and downstream of your target
unigene 

graphication script compute all the different syntenies presente in that unigene clustrer. Thus only show uniques syntenies for every unigene cluster

Highligted in grey appears the ORFs of the unigene cluster, in the central position. Neighbours genes appears in colors, yellowe for Keggs, and green for EGGnogs. In case there are not neighbours genes, these will marked as a red [X]. Functional description are below.
