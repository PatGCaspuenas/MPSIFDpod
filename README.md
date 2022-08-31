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

![utils](http://url/to/img.png)

![data](http://url/to/img.png)



## Acknowledgments

This project has received funding from the European Research Council (ERC) 
under the European Union’s Horizon 2020 research and 
innovation programme (grant agreement No 949085).


