# Simulation of a neuro-epithelium

Vertex model simulation on a toroidal domain, a model for growth and morphogenesis of the developing neural tube.

Python 3 version of this [original code](https://bitbucket.org/Pigueco/vertex_model_python_2.7/src/master/). Study published in

> Guerrero, P., et al. (2019). Neuronal differentiation influences progenitor arrangement in the vertebrate neuroepithelium. *Development*, **146**(23) [https://doi.org/10.1242/dev.176297](https://doi.org/10.1242/dev.176297)

### Requirements

* Python 3.7.9
* numpy
* matplotlib
* numba
* ...


### Installation

It's recommended to setup a conda environment with the above requirements (here called `py37`)
```bash
conda create -n py37 python==3.7.9 numpy matplotlib numba jupyter
conda init $SHELL
conda activate py37
```

Once you setup your python installation and/or activated your conda environment

```bash
git clone https://github.com/berberto/vertex_model.git
cd vertex_model/
pip install -e .
```

### Run

You can test the installation and the functions available with the `vertex_model` from the jupyter notebook in `examples/`:

```bash
jupyter-notebook examples/example_simulation.ipynb
```

#### Main functions

The main functions to look at, in order to see how to set differentiation rates etc, are in `vertex_model/run_select.py`:
- `simulation_with_division`
- `simulation_with_division_clone`
- `simulation_with_division_clone_differentiation`
- `simulation_with_division_clone_differenciation_3stripes`

These are the functions to be called for a simulation run.
The function `run_simulation_INM` is a wrapper for these function, and allows to select which of these to be run, along with other options.