# MPSIFDpod

MPSIFDpod, or Multi-Pulse System Identification of Flow Dynamics with Proper Orthogonal Decomposition, 
is the Matlab implementation of the methodology proposed in the Bachelor Thesis 
"Non-linear system identification for aerodynamic flows" by Patricia García-Caspueñas,  
supervised by Stefano Discetti. This code allows to retrieve a low-order time-resolved velocity 
and pressure fields from incomplete PIV measurements in the temporal domain.

MPSIFDpod uses in input velocity and acceleration snapshot fields to retrieve a reduced-order model 
that serves for inferring a sparse low-order dynamical system of the latter. Once identified, 
the system is integrated with a weighted backward and forward integration process for increased 
accuracy. Estimation of pressure is achieved by enforcing first physical principles, i.e., 
spatial integration of the pressure gradient from Navier Stokes' momentum equation.

## Installation

Installation requires simply that you
[download MPSIFDpod](https://github.com/PatGCaspuenas/MPSIFDpod/archive/refs/heads/main.zip)
and add the base directory (`MPSIFDpod` directory)
to your Matlab path.

### Dependencies

MPSIFDpod has been developed in Matlab R2021b Academic version. Proper working
of the code can be achieved by installing the [cvx package](http://cvxr.com/cvx/).

## Quick usage guide

MPSIFDpod main function `main` is found under the `MPSIFDpod` directory. 
This code serves as an example of the data preparation, training and 
testing processes of the Fluidic Pinball configuration. It can be easily 
generalized for the PIV0018 dataset.

This code requires in input a structure of parameters for some given
functions. These parameters are explained inside each function.
Nevertheless, for a better understanding there exists an interactive
table containing all the required parameters along the process, 
found [here](https://perfect-vibraphone-20c.notion.site/Code-8998ca572e30440f851aedbde6312d18).

MPSIFDpod uses several structures to operate:

* The `TrS` structure contains all parameters corresponding to the training data. That is,
the training dataset, the Proper Orthogonal Decomposition basis and the System Identification
variables that make up for the inferred dynamical system.
* The `SNPM` structure stands for snapshot matrix and includes every general parameter with regards 
to the type of flow and geometry, as well as to the data preparation process. 
* The `TS` structure presents the original time-resolved testing dataset, an undersampled
version of the latter, and the reconstruction models used: Taylor's Hypothesis, cubic
Spline Interpolation and Backward-Forward Integration.

More details on each specific content of each structure is throughly explained inside the dummy-code `main`.

### Code structure

<a href="utils"><img src="https://github.com/PatGCaspuenas/MPSIFDpod/blob/main/docs/utils.png" height="700" ></a>

<a href="datas"><img src="https://github.com/PatGCaspuenas/MPSIFDpod/blob/main/docs/data.png" height="180" ></a>


A directory tree structure of the full repository can be seen above. Some remarks about the repository:

* The `\utils`directory contains all the functions used for the implementation. Arrangement of the functions is seen above. Functions marked with <tspan fill="red">red </tspan> `rgb(255, 0, 0)`
and an asterisk correspond to open source code from external developers. In particular:

    * <font color="red"> Laplacian2D.m </font> and <font color="red"> pIntegrator.m </font> were implemented by Chen, J. and Raiola, M. and Discetti, S. in [Data-driven-pressure-estimation-with-EPOD](https://github.com/erc-nextflow/Data-driven-pressure-estimation-with-EPOD).
    * <font color="red"> poolpolyData.m </font> by S. Brunton, J. Proctor and J. N. Kutz in [MATLAB: SINDy code base](https://faculty.washington.edu/kutz/page26/).
    * <font color="red"> PSD.m </font> by Khmou, Y. in [PSD (Power Spectral Density), and Amplitude Spectrum with adjusted FFT](https://www.mathworks.com/matlabcentral/fileexchange/40002-psd-power-spectral-density-and-amplitude-spectrum-with-adjusted-fft).

* The `\data` directory is for the moment empty in the uploaded repository. It should contain the files depicted above. An 
external link to give access to this data will be uploaded soon. With regards to the raw data necessary for the post-processing and arragement of the datasets:
   * Access to the DNS Fluidic Pinball data was granted by Noack, B.R. and Morzynski, M. in The Fluidic PinBall - a Toolkit for Multiple-Input Multiple-Output Flow Control ({V}ersion 1.0).
    * Access to the planar PIV experiment can be found available in [https://doi.org/10.5281/zenodo.6473075](https://doi.org/10.5281/zenodo.6473075)
    

## Acknowledgments

This project has received funding from the European Research Council (ERC) 
under the European Union’s Horizon 2020 research and 
innovation programme (grant agreement No 949085). We warmly acknowledge M. Morzynski and B. Noack for granting access to the fluidic pinball DNS code, and M. Raiola and J. Chen for providing the 2D wing experimental data and the pressure integration code.


