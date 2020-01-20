# Gene cluster calling
Python scripts for gene cluster calling from the output of OrthoFinder


Software Requirements
---------------------
- OrthoFinder http://www.stevekellylab.com/software/orthofinder


Instructions
------------

You should have one file of protein sequences in fasta format for each species that you plan to include in your analysis. The file names should be short as they will be used during the mirlo analysis to indicate each species in the multiple sequence alignments.

1. run OrthoFinder.

1. run infer_cluster_syntney_from_OG.py. There are 4 input files should be given. 
infer_cluster_syntney_from_OG.py og_group.txt species_map gene_map species_tree > gene_clusters.txt
  og_group.txt is the output file of OrthoFinder: Orthogroups/Orthogroups.txt
  species_map file contains pseudo species with true file name. The sample file located in sample/SpeciesIDs.txt 
  gene_map is the output file of OrthoFinder: WorkingDirectory/SequenceIDs.txt 
  species_tree is the output file of OrthoFinder: Species_Tree/SpeciesTree_rooted.txt


Contact
-------
If you have questions or need help, email me:

Dongqing Yang dongqingyang7@foxmail.com
