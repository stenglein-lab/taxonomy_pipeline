## Create a singularity image that includes R libraries needed for this pipeline

This directory contains singularity definition files to create singularity images containing the R libraries and perl modules needed by scripts in this pipeline.  The R image (r_taxonomy_tools) is based on the [rocker/r-ver docker image](https://hub.docker.com/r/rocker/r-ver) and adds on tidyverse packages and a couple additional packages.  The perl image (perl_taxonomy_tools) is based on an Ubuntu 20.04 docker image and adds some required perl modules.

#### To modify and create the singularity images

Modify the def` files as appropriate

```
# create the singularity image for R scripts
sudo singularity build r_taxonomy_tools.sif r_taxonomy_tools.def

# create the singularity image for perl scripts
sudo singularity build perl_taxonomy_tools.sif perl_taxonomy_tools.def
```

Sign the images:
```
singularity sign r_taxonomy_tools.sif
singularity sign perl_taxonomy_tools.sif
```

If you don't have an existing key you will have to create one using
```
singularity key newpair
```

#### To push this image to the singularity library 

I first had to create a singularity library account and login as described [here](https://sylabs.io/guides/latest/user-guide/cloud_library.html?highlight=push#overview)

After you create a token you will need to do a remote login
```
singularity remote login
```

```
# modify version number as appropriate
singularity push -D "A singularity image containing tidyverse and a few other R packages" r_taxonomy_tools.sif library://stenglein-lab/r_taxonomy_tools/r_taxonomy_tools:1.0.0
singularity push -D "A singularity image containing perl and a few perl modules" perl_taxonomy_tools.sif library://stenglein-lab/perl_taxonomy_tools/perl_taxonomy_tools:1.0.0
```


