# Starfish Discovery

(previously: ilikebigmobileelementsandicannotlie)

## Group members

Jonah Cullen, Callum MacPhillamy, Mauricio Moldes, Alexis Norris, Meghana Ram, Marcus Chan, and Russel Santos

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
5. Wrap R code into Shiny app with `golem` (Jonah)
6. If we discover something in our genomes, follow-up on that :)
