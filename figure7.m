%% Script to recreate Figure 7,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

% --- Create the 6 different stencils

M = 4;
X = [-M:M] -0.5;

N = 0; % Order of differentiation


FD_taylor  = fdweights(0,X,N);
FD_ER0_2   = FD_rel(X, N, pi/2,1e5,1e5,false);
FD_EG0_2   = FD_gv(X, N, pi/2,1e5,1e5,false);
FD_ER0_inf = FD_rel(X, N, pi/2,0,1e5,false);
FD_EG0_inf = FD_gv(X, N, pi/2,0,1e5,false);
FD_hicks0  = FD_hicks( X, M, 6.31 );

ks = linspace(0, 2.5, 10000);
gv_Erf  = @(w,D) abs( 1 - (1-N+1i*w'*X)./(1i*w').^N .* exp(+1i*w'*X) * D' ); % Eq. 15
rel_Erf = @(w,D) abs( 1 - 1./(1i*w').^N .* exp(+1i*w'*X) * D'             ); % Eq. 13
abs_Erf = @(w,D) abs( (1i*w').^N - exp(+1i*w'*X) * D'                     ); % Eq. 10

% Create the plots
figure(1)
subplot(2,1,1)
plot(ks, abs_Erf(ks,FD_taylor) , ...
     ks, abs_Erf(ks,FD_ER0_2)  , ...
     ks, abs_Erf(ks,FD_ER0_inf), ...
     ks, abs_Erf(ks,FD_EG0_2)  , ...
     ks, abs_Erf(ks,FD_EG0_inf), ...
     ks, abs_Erf(ks,FD_hicks0 ))
ylim([0 3e-3])
xlim([0 1.8])
ylabel('|E_{0}|=|E_{r0}|')
title('Total error, relative error')
legend('Taylor coefficients', ...
       '|E_{r0}|_2-optimal' , ...
       '|E_{r0}|_{\infty}-optimal', ...
       '|E_{g0}|_2-optimal' , ...
       '|E_{g0}|_\infty-optimal', ...
       'Hicks (2002)')
grid on


subplot(2,1,2)
plot(ks, gv_Erf(ks,FD_taylor) , ...
     ks, gv_Erf(ks,FD_ER0_2)  , ...
     ks, gv_Erf(ks,FD_ER0_inf), ...
     ks, gv_Erf(ks,FD_EG0_2)  , ...
     ks, gv_Erf(ks,FD_EG0_inf), ...
     ks, gv_Erf(ks,FD_hicks0 ))
ylim([0 3e-3])
xlim([0 1.8])
ylabel('|E_{g0}|')
title('Group-velocity error')
xlabel('k')
grid on
