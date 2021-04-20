%% Script to recreate Figure 4,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

% Run the uniform grid example. After all 3 FD simulations, the created 
% figure '4' will correspond to that in the paper. Note that the figure in
% the paper was processed additionally: a different color scheme was used
% for the image on top:
% https://ch.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap
% and labels for the 3 different cases are added in the second subplot. The
% final text printed in the MATLAB command window will print the average
% simulation time.
example_variable_grid_1D()