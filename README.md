# Starfish Discovery

(previously: ilikebigmobileelementsandicannotlie)

## Group members

Jonah Cullen, Callum MacPhillamy, Mauricio Moldes, Alexis Norris, Meghana Ram, Marcus Chan, and Russel Santos

## Scope of this project

This project aims to improve the usability of [starfish](https://github.com/egluckthaler/starfish). Starfish ([Gluck-Thaler and Vogen. Systematic identification of cargo-mobilizing genetic elements reveals new dimensions of eukaryotic diversity. Nucleic Acids Research 2024.](https://doi.org/10.1093/nar/gkae327)) is a tool for *de novo* giant mobile element discovery/annotation in fungal genomes. 

Starfish currently requires running multiple, separate bash scripts. We aim to create pipelines/workflows for a more accessible, streamlined, and user-friendly implementation of starfish that could expand its use.

Downstream analysis currently involves `R` scripts on [figshare](https://figshare.com/articles/dataset/Supporting_data_for_Systematic_identification_of_cargo-carrying_genetic_elements_reveals_new_dimensions_of_eukaryotic_diversity_/24430447). We aim to improve accessibility of these scripts through implementation in an `R` Shiny app.

Future goals could include the application of Starfish to non-fungal genomes (particularly, mammalian). However, there are anticipated challenges of acquiring the appropriate annotation input files and computational time when moving from small fungal to large mammalian genomes.

## Deliverables

1. Snakemake pipeline for `starfish` (Jonah)  
2. DNAnexus app for `starfish` (Alexis) 
3. Docker workflow for `starfish` (Callum) 
4. R Shiny app for vizualization of `starfish` output (Meghana, Russel)
