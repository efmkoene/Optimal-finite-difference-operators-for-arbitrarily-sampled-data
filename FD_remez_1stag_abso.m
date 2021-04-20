% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

function [D] = FD_remez_1stag_abso(L,grerr)
%FD_REMEZ_1STAG_ABSO   Compute staggered grid optimal first-derivative symmetric coefficients
%   INPUT
%    FD_REMEZ_1STAG_ABSO(L,grerr) computes the coefficients for
%    half-order L, for a first order derivative with absolute cost-function
%  
%   OUTPUT
%    D = L^inf optimal coefficients for half-order scheme
%
%   EXAMPLES
%    FD_remez_1stag_abso(40,1e-5)
%
% Copyright 2019, authors redacted

tic

% --- Number of expected maxima = L (exclude 0)
ks = linspace(1e-7,pi-1e-7,L)';

% --- Dense grid for interpolation
k=linspace(0,pi,L*800);

% --- Matrices D(=LHS) and d(=RHS), eq. (21)
LHS = @(L,E,KS) 2*sin( KS/2 * [1:2:2*L-1] );
RHS = @(L,E,KS) KS + (-1).^([1:L]'+L)*E;

% --- Relative error: eq. (10)
error = @(d,k) 2*sin( k*([1:L+1/2]-1/2) )*d - k;

%/*
% * Run Linf inversion, i.e., the real-valued Remez algorithm 
% */

while 1
    % --- Algorithm step 1: Solve.
    D = LHS(L,grerr,ks) \ RHS(L,grerr,ks);

    % --- Algorithm step 2: Exchange
    intp = error(D,k');

    % --- ........................... Obtain M extrema
    x = abs(intp);
    n = length(x);
    maxind = find(x > [x(1)-1;x(1:n-1)] & x > [x(2:n);x(n)-1]);
    maxind = maxind(1:L);
    errk = intp(maxind);
    newk = k(maxind)';

    % --- Algorithm step 3: Evaluate
    if (max(errk)-abs(grerr))/abs(grerr) < 1e-7
        disp('I have converged.')
        break 
    end 
    ks = newk;

end
toc

figure(1); plot(k, error(D,k') )
hold on
plot(k,k*0+grerr,'k:')
plot(k,k*0-grerr,'k:')
hold off
ylim(4*grerr*[-1 1])

end