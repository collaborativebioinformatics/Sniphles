![logo](./plots/sniphles-logo.png)

# Sniphles

Sniples: Using Sniffles - but phased
- Start from a BAM with phase per read annotated (e.g. from WhatsHap, LongShot)
- Identify phase blocks and split in monophasic (homozygous), diphasic (heterozygous) and unphased blocks
- Loop over haplotype blogs, split BAM using temporary files and make phased bams => good opportunity for massive parallelization
- run sniffles
- filtering
- concatenating the parts
- merging haplotypes back together using SURVIVOR and force-calling Sniffles

BONUS: 
- can handle cram input (Sniffles CAN'T)
- will output a correctly SORTED VCF (Sniffles doesn't care)

## Please cite our work -- here is the ICMJE Standard Citation:

### ...and a link to the DOI:

### You can make a free DOI with zenodo <link>

## Website (if applicable)

## Intro statement

## What's the problem?

## Why should we solve it?

# What is <this software>?

Overview Diagram

# How to use <this software>

# Software Workflow Diagram

![workflow](./plots/SV_Phasing.png)

# File structure diagram
#### _Define paths, variable names, etc_

# Installation options:

### Installing <this software> from Github using Python.

1. `git clone https://github.com/collaborativebioinformatics/Sniphles.git`
2. `cd <repository-folder-path>`
3. `python setup.py install`

# Testing

Python unittests are enabled, run:

```shell script
python -m unittest
```

# Additional Functionality

### DockerFile

<this software> comes with a Dockerfile which can be used to build the Docker image.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd server`
  3. `docker build --rm -t <this software>/<this software> .`
  4. `docker run -t -i <this software>/<this software>`

### Website

There is also a Docker image for hosting the main website. This should only be used for debug purposes.

  1. `git clone https://github.com/NCBI-Hackathons/<this software>.git`
  2. `cd Website`
  3. `docker build --rm -t <this software>/website .`
  4. `docker run -t -i <this software>/website`
