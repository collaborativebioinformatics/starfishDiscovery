# Starfish Discovery

<img src="https://github.com/collaborativebioinformatics/ilikebigmobileelementsandicannotlie/blob/main/894810_cartoon_starfish_with_binoculars_xl-1024-v1-0.png" width="400">

## Group members

Jonah Cullen, Callum MacPhillamy, Mauricio Moldes, Alexis Norris, Meghana Ram, Marcus Chan, and Russel Santos

From 2024 BCM SV Hackathon

<img src="https://github.com/collaborativebioinformatics/starfishDiscovery/blob/main/groupmember_map.png" width="1200">

## Introduction

Mobile genetic elements (also known as transposons or jumping genes) are known to occur across species, from bacteria to humans. Their insertion can have beneficial (e.g., drive evolution), deleterious (e.g., disrupt gene function), or neutral effects to a given genome. Identifying/annotating mobile elements empowers both the understanding of disease and the development of therapeutics. Starfish is a recently developed tool for *de novo* discovery and annotation of giant mobile elements in eukaryotes ([Gluck-Thaler and Vogen. Systematic identification of cargo-mobilizing genetic elements reveals new dimensions of eukaryotic diversity. Nucleic Acids Research 2024.](https://doi.org/10.1093/nar/gkae327)). Specifically, the tool is built for identifying cargo mobile elements (CMEs) in fungal genomes, but the authors claim that _"it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal Starship or a bacterial integrative and conjugative element: a "captain" gene with zero or more "cargo" genes downstream of its 3' end.”_. 

## Methods

### Streamlining starfish pipeline

The current starship analysis (v1.0.0) requires executing individual bash scripts (Figure 1A). To simplify execution, we aimed to create a true pipeline and Docker contatiner (Figure 1B). The Docker container is available to ensure easy execution across operating systems and to ensure reproducibility (see instructions below). A future hackathon can complete the snakemake pipeline (halfway done).

<img src="https://github.com/collaborativebioinformatics/starfishDiscovery/blob/main/starfish_original.png" width="1200">

**Figure 1. Original implementation of starfish.** Excerpted from <https://github.com/egluckthaler/starfish>. It says workflow, but it's individual bash scripts.

### Docker instructions

```bash
# Clone repository
git clone https://github.com/collaborativebioinformatics/ilikebigmobileelementsandicannotlie

# cd into the project repo
cd ilikebigmobileelementsandicannotlie

# Make a new directory for the results to be stored.
mkdir docker_example

# Run the container. 
## Please excuse the personal dockerhub repo. I can move it to a more official one.
docker run -it -v $(pwd)/docker_example:/tmp/starfish_example_output/ callumm93/starfish:v1.0.0

# Inside the container, run the example script.
bash run_starfish_example.sh
```
This should go through the tutorial example and produce some output.

To run `starfish` with new data, use the following:

```bash

# cd to the directory with the data
cd /path/to/your/new/data

docker run -it -v $(pwd):/opt/data --platform linux/amd64 callumm93/starfish:v1.0.0

# The container will put you inside `/opt/starfish/`.
cd ../
ls -1 # You should see data/ here. Inside data/ you should have your new data.
```

### Challenges with creating R Shiny app

We sought to create a Shiny app to visualize/automate the downstream analyses of starfish outputs. This seemed like an easy task given the paper providing code on [figshare](https://figshare.com/articles/dataset/Supporting_data_for_Systematic_identification_of_cargo-carrying_genetic_elements_reveals_new_dimensions_of_eukaryotic_diversity_/24430447?file=44950315) (Figure 2). However, it seems that the figshare was supposed to include directory of additional data (and maybe pre-processing code used to generate it from original starfish output) for each figure, but instead it’s just the pdf of the figure (Figure 2). Without knowing what these 39 files read into the R script look like, it’s hard to match them up to the many output files of starfish. A future hackathon can tackle this issue. For now, we've created a framework for the Shiny app (Figure 3).

<img src="https://github.com/collaborativebioinformatics/starfishDiscovery/blob/main/figshare_missing_data.png" width="1200">

**Figure 2. Missing data.** The `starfishScripts.R` require 39 input files to generate plots. The input file names in the the script are not found in the example run output, suggesting post-processing is required prior to implementing this script. The script references a directory for each figure, suggesting that the authors may have intended to upload a tarball directory for each figure that included the input files and post-processing scripts. However, only the pdf version of each figure is provided in the figshare.

<img src="https://github.com/collaborativebioinformatics/starfishDiscovery/blob/main/shiny_framework.png" width="1200">

**Figure 3. Framework for Shiny app.** Once the `starfishScripts.R` input files are resolved, we have `plotting_for_discovery.Rmd` file with Shiny framework ready to go.

## Conclusions and future directions

Our docker for starfish improves the usability of starfish, as many in our group had difficulties running the program on a Mac. The future completion of the started `snakemake` implementation will streamline the modular `starfish` code to further improve the usability of the tool. Improving the usability will enable the identification and annotation of mobile genetic elements in species beyond fungi. The usability would be dramatically improved through a Shiny app for visualizing the results as the provided code and data from the original study are not reproducible for generating figures. 

Future goals also include the application of `starfish` to non-fungal genomes (particularly, mammalian). However, there are anticipated challenges of acquiring the appropriate annotation input files and computational time when moving from small fungal to large mammalian genomes. An alternative approach could involve using different computational tools to identify transposons in eukaryotic organisms. By employing multiple methods, you can derive a consensus, which may enhance the accuracy and reliability of identifying mobile genetic elements.
