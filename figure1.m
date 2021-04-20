%% Script to recreate that data for Figure 1,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

clear all; close all; clc

N = 1; % Order of differentiation
X = [-6.5:6.5]; % Sample locations

% Create the 5 FD operators. Suppress figures
FD_taylor = fdweights(0,X,N);
FD_ER1_2  = FD_rel(X,N,2.07,1e3,1e3,false);
FD_EG1_2  = FD_gv(X,N,1.715,1e3,1e3,false);
FD_ER1_inf= FD_rel(X,N,2.1,1,1e3,false);
FD_EG1_inf= FD_gv(X,N,1.747,1,1e3,false);

% Plot range (k) and error functions (*_Erf)
ks = linspace(0, 2.5, 10000);
gv_Erf  = @(w,D) abs( 1 - (1-N+1i*w'*X)./(1i*w').^N .* exp(+1i*w'*X) * D' ); % Eq. 15
rel_Erf = @(w,D) abs( 1 - 1./(1i*w').^N .* exp(+1i*w'*X) * D'             ); % Eq. 13
abs_Erf = @(w,D) abs( (1i*w').^N - exp(+1i*w'*X) * D'                     ); % Eq. 10

% Create the plots
figure(1)
subplot(3,1,1)
plot(ks, abs_Erf(ks,FD_taylor) , ...
     ks, abs_Erf(ks,FD_ER1_2)  , ...
     ks, abs_Erf(ks,FD_ER1_inf), ...
     ks, abs_Erf(ks,FD_EG1_2)  , ...
     ks, abs_Erf(ks,FD_EG1_inf), ...
     ks, ks*0+1e-5, 'k--')
ylim([0 3e-5])
xlim([0 2.5])
title('Total error')
legend('Taylor coefficients', ...
       '|E_{r1}|_2-optimal' , ...
       '|E_{r1}|_{\infty}-optimal', ...
       '|E_{g1}|_2-optimal' , ...
       '|E_{g1}|_\infty-optimal')
grid on
ylabel('|E_1|')

subplot(3,1,2)
plot(ks, rel_Erf(ks,FD_taylor) , ...
     ks, rel_Erf(ks,FD_ER1_2)  , ...
     ks, rel_Erf(ks,FD_ER1_inf), ...
     ks, rel_Erf(ks,FD_EG1_2)  , ...
     ks, rel_Erf(ks,FD_EG1_inf), ...
     ks, ks*0+1e-5, 'k--' )
ylim([0 3e-5])
xlim([0 2.5])
title('Relative error')
grid on
ylabel('|E_{r1}|')


subplot(3,1,3)
plot(ks, gv_Erf(ks,FD_taylor) , ...
     ks, gv_Erf(ks,FD_ER1_2)  , ...
     ks, gv_Erf(ks,FD_ER1_inf), ...
     ks, gv_Erf(ks,FD_EG1_2)  , ...
     ks, gv_Erf(ks,FD_EG1_inf), ...
     ks, ks*0+1e-5, 'k--' )
ylim([0 3e-5])
xlim([0 2.5])
title('Group-velocity error')
grid on
ylabel('|E_{g1}|')
xlabel('k')






