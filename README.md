# DWCmod - dropwise condensation modeling

A simple python implementation of thermodynamic models for dropwise condensation. 
Overall heat flux, single drop heat flux and drop size distribution can be calculated for various fluids.
The fluid properties are calculated automatically thanks to [CoolProp](http://www.coolprop.org) \[1\].

At the moment, only one model is implemented. 
Please refer to the [original paper](http://dx.doi.org/10.1115/1.4003742) for a detailed description of the model \[2\]. 


## Getting Started

The easiest way is to use the included jupyter notebook (DWCmod.ipnb) and modify the input values as needed. 
Click here to run the notebook in your browser: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JSablowski/DWCmod/master?filepath=DWCmod.ipynb)

All thermodynamic calculationsc are done in DWC_models.py. DWC_plotting.py contains some basic plotting functions.

The following non-standard python packages are required:
* pylab
* numpy 
* scipy
* coolprop


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


## References

\[1\] Bell, I. H., Wronski, J., Quoilin, S., & Lemort, V. (2014). Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp. *Industrial & Engineering Chemistry Research, 53*(6), 2498–2508. https://doi.org/10.1021/ie4033999

\[2\] Kim, S., & Kim, K. J. (2011). Dropwise Condensation Modeling Suitable for Superhydrophobic Surfaces. *Journal of Heat Transfer, 133*(8), 081502–081502. https://doi.org/10.1115/1.4003742

You can refer to this repository using the following DOI: [![DOI](https://zenodo.org/badge/169960611.svg)](https://zenodo.org/badge/latestdoi/169960611)