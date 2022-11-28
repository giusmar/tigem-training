# Tigem training on bioinformatic tools and technologies

For this trining we can use [Gitpod](https://www.gitpod.io/) as [High-performance computing (HPC)](https://en.wikipedia.org/wiki/High-performance_computing) hosting platofrm. We have structured the [gitpod configuration file](.gitpod.yml) to properly use this github repository on Gitpod, so make sure do not change or delete this file. 

## What you can do

Once you have started the gitpod environment from this repository, the gitpod configuration file automate some steps:
 - Create a docker images for rnaseq using a [docker configuration file](rnaseqdocker). If you want to change or add tools for the simple rnaseq, you have to modify this file.
 - Download a [Rocker](https://rocker-project.org/) images of Rstudio and build a container. Rocker is a community that prebuild images for R with a Rstudio IDE accessible with no-vnc. Rstudio can be accessed with a hosted link by gitpod, found in the Ports page of terminal.
 - Initialize a terminal in the main direcotory.
 
### Slides
[Introduction to R](https://docs.google.com/presentation/d/1OlHEwL4NAAs7pi59HKsmo3aFpR7n5eKnaJLGSw2q0bg/edit?usp=sharing)

### Useful Links
[R](https://www.r-project.org/) [Rstudio](https://posit.co/) [Docker](https://www.docker.com/) [Rocker](https://rocker-project.org/) [Nextflow](https://www.nextflow.io/) [Gitpod](https://www.gitpod.io/) [Nf-core](https://nf-co.re/) [Sequera Nextflow Training](https://training.seqera.io/)
