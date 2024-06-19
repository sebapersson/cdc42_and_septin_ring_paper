# About

GitHub repository for the manuscript **The origin of septin ring size control in budding yeast**.

## Plotting experimental data

The processed experimental data used in the manuscript's plots is located in the `Data` folder. You can find the plotting script in the `Code/process_results/experimental_data` directory. To generate the plots, run any of the R scripts from this directory. The scripts only use standard [Tidyverse](https://www.tidyverse.org/) R packages for plotting.

## Cdc42 and septin ring simulations

### Running simulations

The scripts for the Cdc42 and Septin ring simulations are located in the `Code/cdc42_septin_simulations` directory, while the simulation code for the alternative model can be found in the `Code/cdc42_alt_simulations` directory. The simulation code is in Python, and you can set up the necessary Python environment using the provided Conda YAML file. To install the packages, given a conda installation enter the following command in the terminal:

```bash
conda env create -f cdc42.yml
```

Before running the main simulations, ensure that the pre-simulations specified in the `Run_pre_runs.py` files are completed, these files can also be provided upon request. The models are simulated using a finite-element solver (FEM), with model geometry meshes created in [Gmsh](https://gmsh.info/). A Gmsh installation is not required—since we provide the generated `.msh` file in the simulation directory. When running a simulation the pvd-files and summary statistics (such as cluster area), are stored in the `Results` and `Intermediate` folders, respectively. Given the substantial computational resources required (over 10,000 CPU hours), these simulations are best run on a cluster. We can provide access to the large simulation files (too large for GitHub), upon request.

### Process results

To process results, navigate to the `Code/process_results/cdc42_cluster`, `Code/process_results/cdc42_alt_cluster` or `Code/process_results/septin_ring_size` directories. From the chosen directory run the `plot.R` script. Only standard [Tidyverse](https://www.tidyverse.org/) R packages are used for plotting.

## Particle simulators

### Running simulations

The particle models are simple models designed to investigate various hypotheses about septin ring size regulation. The simulation code is written in Julia, and the models were simulated using Julia version 1.10.1, although the code should be compatible with later Julia versions. To initiate the Julia simulation environment, in this project root directory start Julia, and in the Julia REPL enter:

```julia
] instantiate
```

or alternatively

```julia
import Pkg; Pkg.instantiate
```

The particle simulations can then be run with:

```bash
julia --project=. Code/particle_simulators/run_particle_simulations.jl ARG
```

where `ARG` can be:

* `test_all_diffusion` - run all simulations for the particle model with diffusion.
* `test_all_no_diffusion` - run all simulations for the particle model without diffusion.
* `all_illustrations` - create all illustration plots (illustrate simulators for the MS).

The simulations can take **long** time to run (should be run on a cluster), therefore we can provide intermediate files upon request (as the files are >1GB we cannot put them on GitHub).

### Process results

To process results, navigate to the `Code/process_results/particle_simulators` directory. From this directory run the `plot.R` script. Only standard [Tidyverse](https://www.tidyverse.org/) R packages are used for plotting.