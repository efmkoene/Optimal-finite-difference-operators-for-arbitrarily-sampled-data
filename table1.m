%% Script to recreate that data for Figure 1,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

% The table will be written out to file 'table1.txt'

clear all; close all; clc

N = 0; % Order of differentiation
X = [-2.3:1.7]; % Sample locations

% Create the 10 FD operators. Suppress figures
tab = zeros( 10, length(X) );
tab(1,:)  = fdweights(0,X,N);
tab(2,:)  = FD_hicks(X,3,6.31);
tab(3,:)  = FD_rel(X,N,1,1e3,1e3,false)
tab(4,:)  = FD_gv(X,N,1,1e3,1e3,false)
tab(5,:)  = FD_rel(X,N,1,0,1e3,false)
tab(6,:)  = FD_gv(X,N,1,0,1e3,false)
tab(7,:)  = FD_rel(X,N,2,1e3,1e3,false)
tab(8,:)  = FD_gv(X,N,2,1e3,1e3,false)
tab(9,:)  = FD_rel(X,N,2,0,1e3,false)
tab(10,:) = FD_gv(X,N,2,0,1e3,false)

fileID = fopen('table1.txt','w');
fprintf(fileID,'%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n',tab');
fclose(fileID);