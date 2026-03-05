# Protein Design Pipeline

A minimal RFdiffusion-based pipeline for protein binder generation.

The workflow performs:

RFdiffusion → backbone generation  
ProteinMPNN → sequence design  
AlphaFold2 → structural validation

This repository focuses on the generative stage using RFdiffusion.

## Repository Structure

protein-design-pipeline

scripts/  
    run_rfdiffusion.sh  

configs/  
    design configuration files  

inputs/  
    input structures  

models/  
    RFdiffusion checkpoints  

images/  
    Apptainer container  

outputs/  
    generated designs  

## Running on HPC

Example:

sbatch scripts/run_rfdiffusion.sh configs/headtail.txt

Each SLURM array job generates designs with different random seeds.

## Requirements

RFdiffusion  
SE3Transformer  
PyTorch  
Apptainer / Singularity

## Pipeline Overview

RFdiffusion generates candidate binder backbones constrained by
target interface regions and hotspot residues.

Multiple seeds are sampled in parallel on GPU nodes.

Generated structures can then be passed to downstream sequence
design and structure validation stages.