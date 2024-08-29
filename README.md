# Starfish Discovery

(previously: ilikebigmobileelementsandicannotlie)

## Group members

Jonah Cullen, Callum MacPhillamy, Mauricio Moldes, Alexis Norris, Meghana Ram, and Russel Santos

## Purpose

Extending [starfish](https://github.com/egluckthaler/starfish) for *de novo* giant mobile element discovery/annotation in non-fungi genomes.

## Deliverables

1. Proof of concept that `starfish` works in non-fungi genomes (or doesn’t) and maybe even discover something :)   
2. DNAnexus app for `starfish` 
3. Vizualization of `starfish` output in shiny app (general use -- not specific for our genomes)  

## Steps

![Flowchart v01](https://github.com/collaborativebioinformatics/ilikebigmobileelementsandicannotlie/blob/main/flowchart_2024-08-28-1455.png)


1. Get `starfish` to work (Jonah --> Marcus)  
2. Get data (genome, gene annotation, protein sequence, RepeatMasker annotation) for multiple species (Meghana: human, mouse; Alexis: pig, cow; also Jonah's horse or dog?) 
3. Implement `starfish` code on DNAnexus (create app)
4. Generate R code for loading/viz data, use figshare code as starting point (Russel --> Mauricio --> Alexis) 
5. Wrap R code into Shiny app with `golam` (Jonah)
6. If we discover something in our genomes, follow-up on that :)


## Running `starfish` with Docker

I ran into various issues getting the pipeline to work on a Mac so I went
straight to making a Docker image. 

You should be able to test it with the following steps:

```bash
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