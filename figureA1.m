%% Script to recreate Figure A-1,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

% This script runs the Remez_gv_version.m algorithm, however it contains a
% plotting section on lines 94-107 to plot intermediate results.

clear all; clc; close all;

X = [-5.2:2.8]; % Sample points
N = 1;          % First derivative
K_c=1;          % Optimize up to k=1

% --- Create matrix D: eq. (21)
pm = [ 1i ;
      -1i ];
D =@(w,t) [(1-N+kron( w'*X, pm )) .* exp(kron( w'*X , pm )) , exp( kron( t', pm ) ) ];
% --- Create solution vector b: eq. (21)
b =@(w) kron(w',pm).^N;

% --- Group-velocity error function: eq. (15)
Erf = @(w,D) 1 - (1-N+1i*w'*X)./(1i*w').^N .* exp(+1i*w'*X) * D;

% --- Number of extrema: eq. (22)
M = ceil( (length(X)+1-N)/2 );

% --- Init maxima and phases: step 0 Initialize.
w = linspace( 0, K_c, M );
t = [0:M-1]*pi;

% --- Dense wavenumber vector for interpolation
ks = linspace(0,K_c,2000*length(X));

% --- Initialize error & number of iterations
E = 1e10;
num_its=0;

while 1  
    % --- Create matrix D(=LHS) and d(=RHS)
    LHS = D(w,t); RHS = b(w);
    % --- Cut off the right-most column with phase information
    TMP = LHS(:,end); LHS=LHS(:,1:end-1);
    % --- Do k=0 constraint, eq: (14)
    if w(1)<ks(300);
        LHS(1:2,:) = pm * (X.^N);
        RHS(1:2)   = pm * factorial(N);
    end
    % --- Division by (+- ik)^N, eq. (15)
    LHS = LHS./RHS;
    RHS = RHS./RHS;
    % --- Add additional constraints of eq. (14)
    for j=0:(N-1)
        LHS = [1e2*(X.^j);LHS];       % <--- weigh constraints of eq. 14
        RHS = [0;RHS];
        TMP = [0;TMP];
    end
    LHS = [LHS, TMP];

    % --- Algorithm step 1: Solve. Obtain h, with coefficients D
    C = LHS\RHS;
    C = real( C(1:length(X)) );

    % --- Algorithm step 2: Exchange. First interpolate function
    intp = Erf(ks,C);
    intp(1) = 1 - (X.^N*C)/factorial(N);

    % --- ........................... Obtain M extrema
    x = abs(intp);
    n = length(x);
    maxind = find(x > [x(1)-1;x(1:n-1)] & x > [x(2:n);x(n)-1]);
    maxind( maxind < 10*M ) = [];

    % --- ........................... Include error at k=0 in pool of options, if previous iteration also did!
    maxind = [1;maxind];
            
    % --- ........................... Enforce that the final maximum is at K_c
    if maxind(end)~=length(ks)
        maxind = [maxind;length(ks)];
    end

    % --- ........................... Check if enough maxima exist, otherwise exit the program (happens e.g. with L2 inversion)
    if length(maxind)<M
        disp('Ran out of maxima!')
        break
    end

    % --- ........................... Obtain the 'last' M maxima (can thus also exclude k=0 if it wasn't required!)
    maxind = maxind(end-M+1:end);

    % --- ........................... Obtain associated frequency and angle
    wnew = ks(maxind);
    tnew = angle(intp(maxind))';
    
    % --- Plot
    if num_its < 4
        subplot(2,2,num_its+1)
        plot(ks, abs( Erf(ks,C) ) ,'b' )
        hold on
        plot( w, abs(Erf(w,C)) , 'k.')
        plot( wnew, abs(Erf(wnew,C)), 'rv')
        plot( wnew(1), abs(1 - (X.^N*C)/factorial(N)), 'k.' )
        plot( wnew(1), abs(1 - (X.^N*C)/factorial(N)), 'rv' )
        title(['Iteration ',num2str(num_its+1)])
        grid minor
        ylabel('|E_{g1}|')
        xlabel('k')
    end

    % --- Algorithm step 3: Evaluate: eq. (24).
    Enew = max(abs(intp(10:end)));
    if (abs(E - Enew)/(Enew)) < 1e-7 % Check if error has converged
            break
    end
    if num_its>100
        disp("didn't find it?...")
        break
    end
    E = Enew;

    w = wnew;
    t = tnew;

    num_its=num_its+1;

end

