%% DATASETS OF MPSIFDpod REPOSITORY
%
% Contains fluid data on the flow configurations exploited for the proposed methodology.
% Each configuration is referred with an acronym, such that
%
% 	FP      : Fluidic Pinball                  -> DNS code of M. Morzynski and B. Noack
%	PIV0018 : 2D stalled wing NACA0018 profile -> planar PIV experiment of J. Chena and M.Raiola
%
% There are two main types of files
%
% 	1. X_YYk_ZZk : where X = variable-length string that determines flow configuration, 
%		            YY = number of initial snapshot (*1e3)
%		            ZZ = number of final snapshot (*1e3)
%
% 		       .mat file that contains time-resolved velocity, acceleration and time fields.
% 		       In some cases also pressure and velocity gradients are included. 
%		       Under-sampling must be done afterwards.
%
%		       *Note that the experimental sets are denoted as PIV0018cf, with cf meaning
%			that the whole PIV domain is used as snapshot data and that the latter
%			has been previously filtered with a Saviztky-Golay filter. Acceleration and
%			gradient fields were direcdtly obtained from filtering the velocity.
%
% 	2. X_grid    : where X = variable-length string that determines flow configuration.
%
%		       .mat file containing X and Y grids in non-dimensional form (XD and YD have 
%			SI dimensions, if included). This data is prepared afterwards with the 
%			\+utils\+data_settings\prepareSNPM function
%
%	** Note that NACA0018.txt is a text file with the (x,y) coordinates of the given wing profile. 
% 	   This is only leveraged for visualization purposes
%
