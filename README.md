# The Disproportionate Role of Ocean Topography on the Upwelling of Carbon in the Southern Ocean

### Riley X. Brady, Mathew E. Maltrud, Phillip J. Wolfram, Henri F. Drake, and Nicole S. Lovenduski

Contact: riley.brady@colorado.edu

This repository houses the code used for analysis and visualization in the above manuscript submitted to Geophysical Research Letters. This repository serves to support open-source science by making the analysis transparent and accessible to other scientists. 

The trimmed down Lagrangian particle trajectories can be retrieved at Zenodo. The Zenodo also holds the relavent Eulerian mesh information needed for analysis. 

## How to Use Locally

I have included an `environment.yml` file here. You need to have `miniconda` or `anaconda` installed on your local computer or cluster that you're working on and then run `conda env update -f environment.yml`, which will install all the required packages for this analysis. You can then activate that environment via `conda activate brady-carbonpathways`. 

## Analysis Notebooks

You can find the analysis scripts under the `notebooks` folder. Note that this analysis was run on the Casper cluster at the National Center for Atmospheric Research. The raw particle output is ~50GB and was filtered down to a subset of about 2GB that represents all particle trajectories whose final upwelling across 1000 m occurs S of 45S and outside of the annual Eulerian sea ice margin. 

* `0.00_launch_cluster.ipynb` : Notebook to launch the `dask` cluster, distributed across nodes. The TCP link established here is used to connect other notebooks to the cluster and visualize the progress of parallel workers.
* `1.01_deep_particle_crossing_locations.ipynb` : Post-processes the Lagrangian particle trajectories to find the (x, y) location of the *final* 1000 m crossing for particles. This is used to organize particles into ensembles based on their deep upwelling crossing location.
* `1.02_mixed_layer_particle_crossing_locations.ipynb` : As in the previous notebook, but establishes the (x, y) location of the *first* 200 m crossing for particles following their last 1000 m crossing. This is used to establish ensembles of particles based on their mixed layer crossing location.
* `1.03_1000m_memory_time_origin.ipynb` : Performs the memory time analysis on ensembles of particles to generate post-processed output of the (x, y, z) and tracer concentrations of particles at their statistical origin. I.e., the tracer properties that feed the 1000 m upwelling.
* `1.04_find_mixed_layer_ambient_temperature.ipynb` : Finds the ambient temperature for every particle upon its first crossing into the mixed layer at 200 m. This is used to derive the potential pCO2 values in the paper.
* `1.05_mixed_layer_decorrelation_and_residence_time.ipynb` : Calculates the decorrelation time of DIC and residence time of particles in the mixed layer. This analysis is done based on the 200 m particle ensembles, but evaluates every time a particle crosses into the mixed layer after this first 200 m crossing.
* `2.01_Figure_1.ipynb` : Script and post-processing to create Figure 1 for the paper.
* `2.02_Figure_2.ipynb` : Script and post-processing to create Figure 2 for the paper.
* `2.03_Figure_3.ipynb` : Script and post-processing to create Figure 3 for the paper.
* `2.04_Figure_4.ipynb` : Script and post-processing to create Figure 4 for the paper.
* `3.01_table_results.ipynb` : Script and post-processing to create all of the values reported in Table 1 for the paper.
* `figutils.py` : A script of helpful utility functions to import throughout the notebooks.

**If you are having any trouble viewing these notebooks, paste the URL to the notebook into https://nbviewer.jupyter.org/**.

## Zenodo Data

The following files are provided on the Zenodo server for analysis:

* `bottomDepth.0.5x0.5.nc` : Bathymetry from the Eulerian mesh, regridded onto a 0.5 x 0.5 rectilinear mesh.
* `bottomDepth.nc` : Bathymetry from the Eulerian mesh, on the native unstructured grid.
* `eulerian_sea_ice_climatology.nc` : Monthly climatology of sea ice from the Eulerian mesh, on the native unstructured grid.
* `mesh.nc` : Some mesh descriptors for the unstructured grid, such at latitude, longitude, and cell rea.
* `particle_time_info.nc` : A few helpful formats of time for the time steps on the Lagrangian particles.
* `southern_ocean_deep_upwelling_particles.nc` : Post-processed Lagrangian particle trajectories, including the 19,002 particles that last upwell across 1000 m in the Antarctic Circumpolar Current of the Southern Ocean.

I've also included the post-processed statistical output in `data/postproc` through the Github repository. These post-processed netCDFs can be generated by following analysis included in the notebooks.

## Zarr Format

I provide the Lagrangian particles as a NetCDF file on Zenodo. I tend to use `zarr` these days, which runs a little faster, since it is aware of chunks for `dask`. To convert the NetCDF into a `zarr` format, which I use in the notebooks, do the following:

```python
ds = xr.open_dataset('../data/southern_ocean_deep_upwelling_particles.nc')
ds = ds.chunk({'time': -1, 'nParticles': 'auto'})
ds.to_zarr('../data/southern_ocean_deep_upwelling_particles.zarr', consolidated=True)
```
