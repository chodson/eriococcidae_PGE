We constructed a phylogeny with Eriococcidae species for which we had information on whether heterochromatic bodies (i.e. paternal chrosomes condensed into a heterochromatic ball) were present in male somatic cells. 

## General Information:
- We have two markers available to construct phylogenies with, 18S and CO1. 
- Unfortunately, we don't have sequence data for all the species that we have information on heterochromatic body presence. Also, we have 18S data for more species than we have CO1 data for.
- We made a concantentated phylogeny with all the data we had, which is the phylogeny included in the main text (Figure 1). We also made phylogenies with each of the markers separately which are in the supplementary materials (basically to see if the results we got on the number of losses of heterochromatic bodies was dependent on the marker used). 

## Alignments
- We did alignments in Geneious with the MAFFT plugin. All the information about how we did that is included in the methods of the paper. Below are some additional notes.
- Mask regions of alignment which didn't align properly. Do a strict masking and a less strict masking to see if trees are still congruent. As long as they are use the less strict one.
- Also built phylogenies separately for each gene and see if the tree topology is the same.  

## Notes on making Phylogenies
### Alignments
for CO1 use:
`CO1_Translation_align_MAFFT_Erio_Pnigra_out`
  
for 18S use:
`18S_allseq_17Feb20 for_phylogeny`
  - unedited Mafft phylogeny with all sp. on Lyn's list
  
Concatenated these alignments for the paper:
`Concatenated_align_18sall_co1all_Pnigraout`

### Commands for generating the phylogenies with IQtree
```
qsub -cwd -o logs -e logs -N iqtree -V -pe smp64 4 -b yes 'iqtree -s alignments/18S_allseq_17Feb20_for_phylogeny_alignment.phy -alrt 1000 -bb 1000 -nt AUTO -ntmax 4 -o Parasaissetia_nigra'
qsub -cwd -o logs -e logs -N iqtree -V -pe smp64 4 -b yes 'iqtree -s alignments/Concatenated_align_18sall_co1all_Pnigraout.phy -alrt 1000 -bb 1000 -nt AUTO -ntmax 4 -o Parasaissetia_nigra'
qsub -cwd -o logs -e logs -N iqtree -V -pe smp64 4 -b yes 'iqtree -s alignments/CO1_Translation_align_MAFFT_Erio_Pnigra_out.phy -alrt 1000 -bb 1000 -nt AUTO -ntmax 4 -o Parasaissetia_nigra'
```

**Weird thing** Stictococcus branch is super long. Lyn gave me two stictococcus sequences so I'm adding both to phylogeny to see what it looks like

```
qsub -cwd -o logs -e logs -N iqtree -V -pe smp64 4 -b yes 'iqtree -s alignments/18S_MAFFT_ALL_PNIGRA_STICTO2.phy -alrt 1000 -bb 1000 -nt AUTO -ntmax 4 -o KY927598.1'
```
- looks the same. Read into it, issue is that Stictococcus is not well placed in phylogeny