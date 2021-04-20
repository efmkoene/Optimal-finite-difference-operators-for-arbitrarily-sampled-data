% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

function [C,F,E] = FD_rel(X,N,K_c,extra_maximum,L2_initpoints,display)
%FD_rel   Compute the relative-error optimal coefficients
%   INPUT:
%    FD_rel(X,N,K_c) returns the finite-difference coefficients to
%    get the N-th derivative at x=0 from samples at X, accurate up to a
%    wave-number K_c, with standard number of maxima.
%
%    FD_rel(X,N,K_c,extra_maximum) expects 'extra_maximum' more
%    maxima. For extra_maxima>1, this results in an L2 inversion.
%
%    FD_rel(X,N,K_c,extra_maximum,L2_initpoints) Initializes with
%    L2_initpoints to make an estimate of the extrema and their phases.
%    If you don't want to start with an L2 initialization, choose
%    L2_initpoints < 0, i.e., -1. E.g., FD_rel([-5.2:3.8],1,1,0,-1)
%
%    FD_rel(X,N,K_c,extra_maximum,L2_initpoints,display) allows a 
%    showing of output (default: display=true) which can be turned off (display=false)
%    
%   OUTPUT:
%    D = optimal coefficients,
%    F = Taylor coefficients,
%    E = maximum absolute error.
%
%   EXAMPLES:
%    ASYMMETRIC
%     FD_rel([-3.2:0.8],0,1,0,10)
%     FD_rel([-3.2:0.8]+randn(1,5)*0.2,0,1,0,10)
%     FD_rel([-10.5:10.5]+0.2,0,2.0,0,1e4)
%     X=[-3.2:3.6]+randn(1,7)*0.5;FD_rel(X,3,1,0,500)
%     FD_rel([-10.5:9.5]+0.2,0,1.5,0,1e5)
%    SYMMETRIC
%     FD_rel([-4.5:4.5],1,1.9,1,1e4)
%     FD_rel([-10.5:10.5],2,1.9,0)
%     FD_rel([-10:10],2,1.9,1,1e4)
%
% Copyright 2019, authors redacted

% Default parameters
EM = 0;     % Extra maxima (different from default).
L2 = 1e3;   % Initialize with Least Squares, sampling the k range with 1000 values, eq. (25)
dis=true;   % Show summarizing figure with cost functions

% Optional parameters
switch nargin
    case 4
        EM = extra_maximum;
    case 5
        EM = extra_maximum;
        L2 = L2_initpoints;
    case 6
        EM = extra_maximum;
        L2 = L2_initpoints;
        dis = display;
end

%/*
% * Initialize: set-up all the constants.
% *             set-up the matrices 
% */

% --- Create matrix D: eq. (21)
pm = [ 1i ;
      -1i ];
D =@(w,t) [exp( kron( w'*X, pm ) ) , exp( kron( t', pm ) ) ];
% --- Create solution vector b: eq. (21)
b =@(w) kron(w',pm).^N;

% --- Relative error function: eq. (12)
Erf = @(w,D) 1 - 1./(1i*w').^N .* exp(+1i*w'*X) * D;

% --- Number of extrema: eq. (22)
M = ceil( (length(X)+1-N)/2 );
M = M+EM; % Add additional tested extremal wavenumbers

% --- Init maxima and phases: step 0 Initialize.
w = linspace( 0, K_c, M );
t = [0:M-1]*pi;

% --- Dense wavenumber vector for interpolation
ks = linspace(0,K_c,2000*length(X));

%/*
% * Run L2 inversion, used to initialize the maxima & phases
% */

if EM<=1 && L2>0
    % --- Obtain initial estimate of maxima and phases w and t
    disp('Obtaining initial estimate using L2 inversion')
    C = FD_rel(X,N,K_c,L2,[],false)'; % Run L2-inversion
    intp = Erf(ks,C);
    intp(1) = 1 - (X.^N*C)/factorial(N);
    % --- Obtain M extrema
    x = abs(intp);
    n = length(x);
    maxind = find(x > [x(1)-1;x(1:n-1)] & x > [x(2:n);x(n)-1]);
    if maxind(end)~=length(ks)
        maxind = [maxind;length(ks)];
    end
    maxind = maxind(end-M+1:end); % Obtain 'last' values...
    % --- Set initial extrema and phases
    w = ks(maxind(end-M+1:end));
    t = angle(intp(maxind(end-M+1:end)))';
    %--------------------------------------------------------------- //end
end


%/*
% * Run Linf inversion, i.e., the Complex Remez algorithm 
% */

% --- Initialize error & number of iterations
E = 1e10;
num_its=0;

tic
while 1  
    % --- Create matrix D(=LHS) and d(=RHS)
    LHS = D(w,t); RHS = b(w);
    % --- Cut off the right-most column with phase information
    TMP = LHS(:,end); LHS=LHS(:,1:end-1);
    % --- Do k=0 constraint, eq: (14)
    if w(1)<ks(300)
        LHS(1:2,:) = pm * (X.^N);
        RHS(1:2)   = pm * factorial(N);
    end
    % --- Division by (+- ik)^N, eq. (15)
    LHS = LHS./RHS;
    RHS = RHS./RHS;
    % --- Add additional constraints of eq. (14)
    for j=0:(N-1)
        LHS = [5e1*(X.^j);LHS];       % <--- weigh constraints of eq. 14
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
toc
disp(['~~~~~~~num_its=',num2str(num_its),' with error=',num2str(max(abs(intp))), '. Expected maxima=',num2str(M)])

% --- Approx target (Taylor coefficients)
F = fdweights(0,X,N);

%/*
% * The program is strictly done now, we'll just plot some things
% */
if dis==true
    Erf = @(w,D) 1 - (1-N+1i*w'*X)./(1i*w').^N .* exp(+1i*w'*X) * D; % Eq. 15
    rel_Erf = @(w,D) 1 - 1./(1i*w').^N .* exp(+1i*w'*X) * D; % Eq. 13
    abs_Erf = @(w,D) (1i*w').^N - exp(+1i*w'*X) * D;         % Eq. 10
    
    figure(1);
    clf
    subplot(2,2,1)
    plot(ks, real(Erf(ks,C)),'r', ks, imag(Erf(ks,C)),'b')
    hold on
    plot( w, real( Erf(w,C) .* exp(+i*t)' ) , '*')
    plot(ks, abs( Erf(ks,C) ) ,'g' )
    plot(ks, real(Erf(ks,F')),'r:', ks, imag(Erf(ks,F')),'b:')
    hold off
    grid on
    title('Group-velocity error')
    
    subplot(2,2,2)
    plot(ks, real(rel_Erf(ks,C)),'r', ks, imag(rel_Erf(ks,C)),'b')
    hold on
    plot( w, real( rel_Erf(w,C) .* exp(+i*t)' ) , '*')
    plot(ks, abs( rel_Erf(ks,C) ) ,'g' )
    plot(ks, real(rel_Erf(ks,F')),'r:', ks, imag(rel_Erf(ks,F')),'b:')
    hold off
    grid on
    title('Relative error')
    
    subplot(2,2,4)
    plot(ks, real(abs_Erf(ks,C)),'r', ks, imag(abs_Erf(ks,C)),'b')
    hold on
    plot( w, real( abs_Erf(w,C) .* exp(+i*t)' ) , '*')
    plot(ks, abs( abs_Erf(ks,C) ) ,'g' )
    plot(ks, real(abs_Erf(ks,F')),'r:', ks, imag(abs_Erf(ks,F')),'b:')
    hold off
    grid on
    title('Absolute error')
    
    subplot(2,2,3)
    stem(X,[F' C])
    legend('Taylor','Optimal')
    title('Space-domain FD coefficients')
end

C=C';
end