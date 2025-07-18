# RNA-seq R docker

This directory contains the recipe to build the docker container for R analysis of RNA-seq data. This uses the bioconductor docker as a base so you can port to using `Rstudio` in the browser.

## Downloading the image
To download, the image, you can simply run

```bash
docker pull kwellswrasman/rna_seq_r:v2
```

## Converting to singularity

To convert this image to singularity, simply run
```bash
singularity pull --name rna_seq_r_v2.sif docker://kwellswrasman/rna_seq_r:v2
```

## Using R studio

To run `Rstudio` using this image, in a terminal run

```bash
docker run \
	-e PASSWORD=bioc \
	-p 8787:8787 \
	kwellswrasman/rna_seq_r:v2
```

If you navagate to http://localhost:8787 on you web browswer, you will be able to log into a Rstudio session using `rstudio` as your username and whatever you set above as the password, in this case `bioc`.

With the above command, you won't have access to your system, but adding a mount line will fix it.

```bash
docker run \
	-e PASSWORD=bioc \
	-p 8787:8787 \
	--mount type=bind,source="$(pwd)",target=/home/rstudio/rnaseq \
	kwellswrasman/rna_seq_r:v2
```

More information on the base container is [here](https://bioconductor.org/help/docker/)

Packages installed in this container can be found in the `R_dependencies` file. The `renv.lock` file will provide all packages and versions.

## R studio with singularity
To run R studio with singularity on a slurm server, use the helper `launch_rstusio.sh` script. The log file will include instruction for how to run rstusio from within the singularity image in a running job. This requires ssh access to the server.

## Version control with `renv`
This image was built using version control with `renv`. 

### Initial build
To build the container initially, I added any desired packages to `R_dependencies` and then I followed the following steps

1. I removed the following lines from the docker file

```
COPY renv.lock renv.lock
RUN R -e "renv::restore()"
```

2. I then built the package without any R packages
```bash
docker build rna_seq_r:v2 ./
```

3. I next started an interactive shell

```bash
docker run -it --mount type=bind,source="$(pwd)",target=/home/rstudio/r_docker rna_seq_r:v2 sh
```

4. In this I started R and used `renv` to install packages. Any non-cran packages need to be installed manually using the full github path or `bioc::` for bioconductor packages. I then copy the lock file into the r_docker directory.

```bash
R

> renv::init()
> renv::hydrate()
> renv::install(c("github_user/github_package", "bioc::bioconductor_package"))
> renv::snapshot
> q()

cp renv.lock r_docker
exit
```

5. Now that I have the lock file, I add back in two lines from the docker file
```
COPY renv.lock renv.lock
RUN R -e "renv::restore()"
```

6. Rebuild the image
```bash
docker build atac_chip_r:v1 ./
```

### Adding a package

1. Add you new package to `R_dependencies`

2. Start an interactive shell in the docker container (this assumes it has been downloaded from dockerhub, see above)

```bash
docker run -it --mount type=bind,source="$(pwd)",target=/home/rstudio/r_docker rna_seq_r:v2 sh
```

3. Start R and install your new packages with `renv`. Any non-cran packages need to be installed manually using the full github path or `bioc::` for bioconductor packages. I then copy the lock file into the r_docker directory.

```bash
R

> renv::init()
> renv::hydrate()
> renv::install(c("github_user/github_package", "bioc::bioconductor_package"))
> renv::snapshot
> q()

cp renv.lock r_docker
exit
```

4. Rebuild the image
```bash
docker build rna_seq_r:v2 ./
```