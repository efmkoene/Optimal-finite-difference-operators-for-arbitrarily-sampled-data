%% Script to recreate that data for Figure 1,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

% The table will be written out to file 'table3.txt'

clear all; close all; clc

N = 1; % Order of differentiation
X = [-18.5:18.5]; % Samples
L = length(X)/2;  % Operator half-order

% Create the 8 FD operators.
tab = zeros( 8, 2*L );
tab(1,:)       = fdweights(0,X,N);
tab(2,:)       = FD_tot(X, N, 2.882,1e4);
tab(3,L+1:end) = FD_remez_1stag_abso(L,1e-4)';
tab(4,:)       = FD_rel(X, N, 2.929,1e4);
tab(5,L+1:end) = FD_remez_1stag_relative(L,1e-4)';
tab(6,:)       = FD_gv(X, N, 2.692,1e5);
tab(7,L+1:end) = FD_remez_1stag_groupvel(L,1e-4)';
tab(8,L+1:end) = [0.1271216E+1 -0.1394578E+0 0.4893614E-1 -0.2402039E-1 0.1379334E-1 -0.8643880E-2 0.5709894E-2 -0.3896436E-2 0.2711002E-2 -0.1905154E-2 0.1342289E-2 -0.9420260E-3 0.6543898E-3 -0.4468455E-3 0.2973663E-3 -0.1905168E-3 0.1150882E-3 -0.6229052E-4 0.1996929E-4];

% Truncate to half of the operator
tab = tab(:,L+1:end);

fileID = fopen('table3.txt','w');
fprintf(fileID,'Taylor weights\t|E_1|_2-optimal\t|E_1|_inf-optimal\tE_r1|_2-optimal\t|E_r1|_inf-optimal\t|E_g1|_2-optimal\t|E_g1|_inf-optimal\tLiu (2014)\n')
fprintf(fileID,'%1.15f\t%1.15f\t%1.15f\t%1.15f\t%1.15f\t%1.15f\t%1.15f\t%1.15f\n',tab);
fclose(fileID);