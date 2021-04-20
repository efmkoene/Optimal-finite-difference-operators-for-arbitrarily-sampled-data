%% Script to recreate that data for Figure 2,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

function [] = example_uniform_grid_1D()

dx = 10;
N  = 301;
X  = 0:dx:(N-1)*dx;
dt = 0.001;
T  = 0:dt:0.8;
nt = length(T);
c  = 2000;

% --- pre-Compute Forward Time Dispersion Transform, cf. Koene et al. (2018); Eliminating time dispersion from seismic wave modeling; GJI
A= exp( -2i * sin([0:nt-1]*pi/(2*nt))' * [0:nt-1] );
B = ifft( A, 2*nt, 'symmetric');
FTDT = B(1:nt,1:nt);

% --- pre-Compute Inverse Time Dispersion Transform, cf. Koene et al. (2018); Eliminating time dispersion from seismic wave modeling; GJI
A1=exp( -2i * sin([0:nt-1]*pi/(2*nt))' * [0:2*nt-1] );
v  = cos([0:nt-1]*pi/(2*nt));
A2 = v' .* A1;
B  = ifft( A2, 2*nt, 'symmetric' );
ITDT  = B(1:nt,1:nt)';

% --- Wavelet
fc = 40;          % Choose 25 Hz or 50 Hz (or anything else!)
ricker = @(fm,t) (1-2*pi^2*fm^2*t.^2) .* exp(-pi^2*fm^2.*t.^2);
% Input
fs0= ricker(fc, T-0.1 );
% Transformed
fs = FTDT*fs0(:);   

%% Do the modeling for a fixed error of 1e-4
% --- Create CENTERED, STAGGERED finite-difference stencil (order = L)
L = 19;
Xur= [-L:L-1];

for var=1:8
    switch var
        case 1
            % --- Case 1: Taylor
            a = fdweights(0,Xur+0.5,1);
        case 2
            % --- Case 2: Total L_2
            a= FD_tot(Xur+0.5, 1, 2.882,1e4);
        case 3
            % --- Case 3: Total L_inf
            a = FD_remez_1stag_abso(L,1e-4)'; a=[-a(end:-1:1) a];
        case 4
            % --- Case 4: Relative L_2
            a = FD_rel(Xur+0.5, 1, 2.929,1e4);
        case 5
            % --- Case 5: Relative L_inf
            a = FD_remez_1stag_relative(L,1e-4)'; a=[-a(end:-1:1) a];
        case 6
            % --- Case 6: Groupvel L_2
            a = FD_gv(Xur+0.5, 1, 2.692,1e4);
        case 7
            % --- Case 7: Groupvel L_inf
            a = FD_remez_1stag_groupvel(L,1e-4)'; a=[-a(end:-1:1) a];
        case 8
            % --- Case 8: Liu (2014), Optimal staggered-grid ﬁnite-diﬀerence schemes based
            % on least-squares for wave equation modelling: Geophysical Journal International, 197, 1033–1047)
            a = [0.1271216E+1 -0.1394578E+0 0.4893614E-1 -0.2402039E-1 0.1379334E-1 -0.8643880E-2 0.5709894E-2 -0.3896436E-2 0.2711002E-2 -0.1905154E-2 0.1342289E-2 -0.9420260E-3 0.6543898E-3 -0.4468455E-3 0.2973663E-3 -0.1905168E-3 0.1150882E-3 -0.6229052E-4 0.1996929E-4]; a=[-a(end:-1:1) a];
    end
% --- 'Setup' physical wave quantities
p = zeros(N+2*L,1);
v = zeros(N+2*L,1);

% --- Updateable portion
Xu = [L:N+L-1];

% --- Setup recording
rec = zeros(nt,1);

% Create wave simulation
figure(1)
clf
tic
for t = 1:nt   
    % Record P
    rec(t) = p(150);
    
    % Inject p
    p(75) = p(75) + dt*c*2/dx*fs(t); % B
    
    % Update P & V
    p(Xu+1) = p(Xu+1) - dt*c^2*v(Xur+Xu'+1) * a'/dx;
    v(Xu  ) = v(Xu  ) - dt*    p(Xur+Xu'+1) * a'/dx;
    
    % Apply boundary conditions
    %     --- Free surface at left boundary
    p(L) = 0;
    p(1:L-1) = -p(2*L-1:-1:L+1);
    v(1:L)   =  v(2*L:-1:L+1);
    %     --- Free surface at right boundary
    p(N+L) = 0;
    p(N+L+1:N+2*L) = -p(N+L-1:-1:N);
    v(N+L+1:N+2*L) =  v(N+1:-1:N+1);
    
    plot( X, p(Xu))
    title(['case ',num2str(var),'/15'])
    hold on
    plot(X(150-L-1),p(150),'o')
    hold off
    drawnow
end
computational_time(var) = toc;
recs(:,var) = rec;

end

%% --- case 9: Staggered pseudospectral solution
var=9;
N2 = length(p);
dF = 1/(N2*dx);
kx = 2*pi*dF*(0:N2-1)'; % Wavenumber k for pseudospectral computations

% --- 'Setup' physical wave quantities
p = zeros(N2,1);
v = zeros(N2,1);

figure(1)
tic
for t = 1:nt   
    % Record P
    rec(t) = p(150);
    
    % Inject p
    p(75) = p(75) + dt*2*c/dx*fs(t);
    
    % Update P
    tmp2 = ifft( 1i*kx .* exp(-kx*1i*dx/2) .* fft(v),'symmetric'); % vx derivative
    p(Xu)   = p(Xu)   - dt*c^2*tmp2(Xu);
    % Update v
    tmp4 = ifft( 1i*kx .* exp(-kx*1i*dx/2) .* fft(p),'symmetric'); % px derivative
    v(Xu-1) = v(Xu-1)- dt*tmp4(Xu);
    
    plot( X, p(Xu))
    title(['case ',num2str(var),'/15'])
    hold on
    plot(X(150-L-1),p(150),'o')
    hold off
    drawnow
end
computational_time(var) = toc;

recs(:,var) = rec;

for i=1:var
    rectf(:,i) = ITDT * (recs(:,i));
end

figure(2)
subplot(2,2,2)
plot( T, rectf(:,[1:9])-[11:19]*0.45 )
xlim([0.3 0.65])
ylim([-9.5 -1])
legend('Taylor coefficients','|E_1|_2-optimal','|E_1|_\infty-optimal',...
       '|E_{r1}|_2-optimal','|E_{r1}|_\infty-optimal','|E_{g1}|_2','|E_{g1}|_\infty', 'Liu','Pseudospectral','Location','NorthEast')
title('|E|=10^{-4}')
view([90 90])

subplot(2,2,4)
plot( T, 20*(rectf(:,[1:8])-rectf(:,9))-[11:18]*0.45 )
xlim([0.3 0.65])
ylim([-9.5 -1])
% legend('Taylor coefficients','Total L2','Total Linf',...
%        'Relative L2','Relative Linf','Groupvel L2','Groupvel Linf', 'Relative Liu','Location','NorthWest')
title('Error compared to pseudospectral solution \times 20')
view([90 90])

   
%% Case 2: Up to Kc=2.9
for var=10:15
    switch var
        case 10
            % --- case 10: Total L_2
            a= FD_tot(Xur+0.5, 1, 2.9,1e4);
        case 11
            % --- case 11: Total L_inf
            a = FD_remez_1stag_abso(L,9.4e-5)'; a=[-a(end:-1:1) a];
        case 12
            % --- case 12: Relative L_2
            a = FD_rel(Xur+0.5, 1, 2.9,1e4);
        case 13
            % --- case 13: Relative L_inf
            a = FD_remez_1stag_relative(L,4.1e-5)'; a=[-a(end:-1:1) a];
        case 14
            % --- case 14: Groupvel L_2
            a = FD_gv(Xur+0.5, 1, 2.9,1e4);
        case 15
            % --- case 15: Groupvel L_inf
            a = FD_remez_1stag_groupvel(L,4e-3)'; a=[-a(end:-1:1) a];
    end

    % --- 'Setup' physical wave quantities
p = zeros(N+2*L,1);
v = zeros(N+2*L,1);

% --- Updateable portion
Xu = [L:N+L-1];

% --- Setup recording
rec = zeros(nt,1);

% Create wave simulation
figure(1)
clf
tic
for t = 1:nt   
    % Record P
    rec(t) = p(150);
    
    % Inject p
    p(75) = p(75) + dt*c*2/dx*fs(t); % B
    
    % Update P & V
    p(Xu+1) = p(Xu+1) - dt*c^2*v(Xur+Xu'+1) * a'/dx;
    v(Xu  ) = v(Xu  ) - dt*    p(Xur+Xu'+1) * a'/dx;
    
    % Apply boundary conditions
    %     --- Free surface at left boundary
    p(L) = 0;
    p(1:L-1) = -p(2*L-1:-1:L+1);
    v(1:L)   =  v(2*L:-1:L+1);
    %     --- Free surface at right boundary
    p(N+L) = 0;
    p(N+L+1:N+2*L) = -p(N+L-1:-1:N);
    v(N+L+1:N+2*L) =  v(N+1:-1:N+1);
    
    plot( X, p(Xu))
    title(['case ',num2str(var),'/15'])
    hold on
    plot(X(150-L-1),p(150),'o')
    hold off
    drawnow
end
computational_time(var) = toc;


recs(:,var) = rec;
end

% --- Postprocess
for i=10:15
    rectf(:,i) = ITDT * (recs(:,i));
end

figure(2)
subplot(2,2,1)
plot( T, rectf(:,[1 10:15 8 9])-[11:19]*0.45 )
xlim([0.3 0.65])
ylim([-9.5 -1])
legend('Taylor coefficients','|E_1|_2-optimal','|E_1|_\infty-optimal',...
       '|E_{r1}|_2-optimal','|E_{r1}|_\infty-optimal','|E_{g1}|_2','|E_{g1}|_\infty', 'Liu','Pseudospectral','Location','NorthEast')
title('K_c=2.9')
view([90 90])

subplot(2,2,3)
plot( T, 20*(rectf(:,[1 10:15 8])-rectf(:,9))-[11:18]*0.45 )
xlim([0.3 0.65])
ylim([-9.5 -1])
title('Error compared to pseudospectral solution \times 20')
view([90 90])

disp(['The average simulation time = ',num2str(mean(computational_time)),'seconds'])

end