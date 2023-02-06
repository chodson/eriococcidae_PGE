We wanted to determine in _Cystococcus_ species, whether male parents transmit both alleles they carry to offspring, or only one allele (as you would expect if males exhibit PGE since they would only transmit the maternally derived chromosomes/alleles to offspring). We looked at this in two species:
- _Cystococcus echiniformis_: which we found did not carry heterochromatic bodies and so may lack PGE transmission
- _Cystococcus campanidorsalis_: which we found does carry heterochromatic bodies and so we expect to exhibit PGE transmission similar to other Eriococcidae species.

## General notes:
- We collected galls for this analysis. In Cystococcus species, females form galls which she then lives in and lays all her offspring inside, producing male offspring first followed by female offspring. We used galls containing the mother and male offspring for the analysis. This means that we did not know the genotype of the father of the offspring, and had to infer this from the data. We also had to infer whether the female mated once or multiple times.
- We used microsatellite markers for this analysis. We analysed 9 loci for _C. echiniformis_ and 8 loci for _C. campanidorsalis_.


## Primer design:
we designed primers for this analysis by identifying microsatellite loci from low coverage shotgun sequencing data. We used QDD <https://people.imbe.fr/~emeglecz/qdd> to design primers, after we trimmed the reads and assembled the reads with clc.

## Analysis of allele inheritance patterns:
Analysed in R with script `inheritance_analysis.R`