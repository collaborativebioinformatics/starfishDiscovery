# Exploring the fungi data

1. Download the data from https://figshare.com/articles/dataset/Supporting_data_for_Giant_Starship_elements_mobilize_accessory_genes_in_diverse_fungal_genomes_/17185880?file=31767452
2. `gunzip supporting_data.tar.gz `
3. ` tar -xvf supporting_data.tar `
4. copy and rename the assemblies > assembly and gene_gff > gff3 for the bash file
`cp -r assemblies /mnt/c/Users/.../ilikebigmobileelementsandicannotlie/examples/test/assembly`
`cp -r gene_gff /mnt/c/Users/.../ilikebigmobileelementsandicannotlie/examples/test/gff3`
<!-- `cp -r assemblies assembly` -->
<!-- `cp -r gene_gff gff3` -->
<!-- `GFF_PATH="$(pwd)/gff"` -->
<!-- `ASM_PATH="$(pwd)/assembly"` -->

4. Run docker with the file

    - Open Docker Desktop Settings:

    - Right-click on the Docker icon in the system tray and select "Settings."
        WSL 2 Integration:

        In the settings menu, navigate to the “Resources” tab, then select “WSL Integration.”
        Ensure that your WSL 2 distributions (such as Ubuntu) are enabled for Docker integration.
        Apply and Restart:

    - Apply the changes and restart Docker Desktop if necessary.


```
# Clone repository
git clone https://github.com/collaborativebioinformatics/ilikebigmobileelementsandicannotlie

# cd into the project repo
cd ilikebigmobileelementsandicannotlie

# Make a new directory for the results to be stored.
mkdir docker_example

# Run the container. 
## Please excuse the personal dockerhub repo. I can move it to a more official one.
sudo docker run -it -v $(pwd)/docker_example:/tmp/starfish_example_output/ callumm93/starfish:v1.0.0

# Inside the container, run the example script.
bash run_starfish_example.sh
This should go through the tutorial example and produce some output.

To run starfish with new data, use the following:

# cd to the directory with the data
cd /path/to/your/new/data

docker run -it -v $(pwd):/opt/data --platform linux/amd64 callumm93/starfish:v1.0.0

# The container will put you inside `/opt/starfish/`.
cd ../
ls -1 # You should see data/ here. Inside data/ you should have your new data.

# To export data out 
## Open a new terminal 
## Find the Container ID:
docker ps
## Use the docker cp command to copy files from the container to your local machine
docker cp <container_id>:/path/in/container /path/on/host
e.g.  docker cp 45a63f64bd8e:"/opt/starfish/test" "C:\Users\TechD\...\ilikebigmobileelementsandicannotlie\examples\test_fungi"

# Upload to dna nexus
dx upload -r "C:/Users/.../ilikebigmobileelementsandicannotlie/examples/test_fungi"
```
