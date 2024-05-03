# About

## Plotting experimental data

The processed experimental data used in the manuscript's plots is located in the `Data` folder. You can find the plotting script in the `Code/process_results/experimental_data` directory. To generate the plots, run any of the R scripts from this directory. The scripts only use standard [Tidyverse](https://www.tidyverse.org/) R packages for plotting.

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