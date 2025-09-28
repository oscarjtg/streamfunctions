# streamfunctions
Python module for calculating and visualising stream functions from velocity components on a staggered C-grid

## Installation instructions

Clone the repository, and (preferably in a virtual environment) run
```
pip install <path-to-directory>
```
where `<path-to-directory>` should be replaced by the path to the `streamfunctions` directory that you have just cloned (where this project's `pyproject.toml` file is located).

Alternatively, just clone the repository and add it to your Python path, 
or simply copy the module file `streamfunctions.py` to the directory containing the scripts in which you would like to import the module.

## Quick use

Run
```python
import streamfunctions as stfs
```

Load in some two-dimensional flow data computed by a finite volume method on a staggered grid and compute the streamfunction using 
```python
stfs.streamfunction_direct_integration(u, v, xC, xF, yC, yF, psi_bot, psi_left)
```
where 
* u is a two-dimensional NumPy array of shape `(m, n-1)` containing the x-component of the velocity field at the vertical faces of the finite volume grid cells,
* v is a two-dimensional NumPy array of shape `(m-1, n)` containing the y-component of the velocity field at the horizontal faces of the finite volume grid cells,
* xC is a one-dimensional NumPy array of length `(m-1)` containing the x-coordinates of the cell centres,
* xF is a one-dimensional NumPy array of length `(m)` containing the x-coordinates of the cell centres,
* yC is a one-dimensional NumPy array of length `(n-1)` containing the y-coordinates of the cell centres,
* yF is a one-dimensional NumPy array of length `(n)` containing the y-coordinates of the cell centres,
* psi_bot is a one-dimensional NumPy array of length `(m)` containing the bottom boundary values of the stream function, and
* psi_bot is a one-dimensional NumPy array of length `(n)` containing the left boundary values of the stream function.

## TODO
* Improve the calculation by computing the stream function by first computing the vorticity and then solving Poisson's equation for the stream function with the vorticity as the source term.
* Add support for Neumann boundary conditions.
* Quantitatively compare the direct integration and Poisson calculations, for both idealised problems with known analytic solutions and for complex fields from ocean model run output.
