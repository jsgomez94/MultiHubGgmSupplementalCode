# HubEstimationCodeSubmission
This repository contains all the proof-read and commented code to recreate simulations, figures and tables of the paper "Hub Detection in Large Gaussian Graphical Models."

Here, we simply describe the overall structure of the repository. To fully reproduce the results, further instructions are necessary. For this end, we provide further README files when necessary. The structure is the following:

- <code>./figures_4/</code>: directory containing the code for generating all figures that do not require extensive simulations or real data. Each figure has its own directory, and contains individual README files describing how to replicate each of the figures. The code for all figures can be replicated in a personal computer, except Figure S3. For Figure S3, it is necessary to run simulations on a linux cluster with SLURM scheduling system.

- <code>./realdata_4/</code>: directory containing the code and data for performing the real data analysis we describe in Section 6 of our main paper. We provide a README file providing further instructions on how to replicate the real data analysis.

- <code>./simulations_4/</code>: directory containing the code for reproducing the simulation results of Section 5 in our main paper, and Section S3.4 of our Supplementary Materials. In order to replicate all simulation experiments and plots, you must run the simulation code on a linux cluster with SLURM scheduling system. We provide further description of how to reproduce the results in the README of this directory.