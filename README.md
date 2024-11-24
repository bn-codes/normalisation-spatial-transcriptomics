# Normalisation algorithms for spatially resolved transcriptomic data
This repository contains the source code of my bachelor thesis in Bioinformatics at the Faculty of Science, Charles University in Prague.

The thesis is freely available [here](https://dspace.cuni.cz/bitstream/handle/20.500.11956/181526/130358464.pdf).

This program aims to allow the comparison of normalisation algorithms for spatially resolved transcriptomic data. Its output is a visualisation of the normalised data through graphs, and the interpretation must be done manually.

## Basic Structure

The program has nine parts that must be run in a specified order. Seven parts of the program are carried out in R, while the remaining two must be executed in Python. Data transfer between these two languages is provided within the program.

All code is provided in specified order in the folder `Code/`.

 1. Libraries and Environment [R] loads libraries and sets up an environment
 2. Scanpy Normalisation [Python] executes the Scanpy normalisation algorithm
 3. stLearn Normalisation [Python] executes the stLearn normalisation algorithm
 4. Specify Paths [R] determines the remaining necessary paths
 5. Load and Normalise [R] executes the Giotto, Seurat and scran normalisation algorithms; loads the Scanpy and stLearn normalised data obtained in Python; loads clusters of spots
 6. Spot and Gene Groups [R] divides genes and spots into their respective groups
 7. Supplementary Data [R] creates a visualisation of supplementary data
 8. Visualisation Part 1 [R] creates a visualisation of (not)-scaled normalised data based on logarithmic spot UMI count for raw data and all normalisation algorithms
 9. Visualisation Part 2 [R] creates a visualisation of the variance contribution of spots within gene groups

## Input Data Representation

The main input of the program is the RDS input file, which contains un-normalised information about the transcriptomic data obtained from the sample.

It is also necessary to have access to the usual 'outs' folder, the result of 10X Genomics Space Ranger, and spot clusters (the clusters used in this project are provided in `Clusters/`).

## Output Data Representation

The result of the program is a number of graphs:

 1. Graph showing the number of genes in gene groups
 2. Graph showing the number of spots in spot groups
 3. Graph showing the variance contribution of spots within the sample
 4. Graphs showing a visualisation of (not)-scaled normalised data based on logarithmic spot UMI count for raw data and all normalisation algorithms
 5. Graphs showing a visualisation of variance contribution of spots within gene groups for raw data and all normalisation algorithms

The plots created with the above-specified clusters are located in the folder `Output/`.

## Specification of Paths

Several file paths must be specified for the program to run correctly. This part of the code is always indicated by a comment that starts with "specify path", and occurs in program sections 1, 2, 3 and 4.
