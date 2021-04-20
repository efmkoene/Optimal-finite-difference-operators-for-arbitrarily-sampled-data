%% 1D FINITE-DIFFERENCE PROPAGATION WITH ARBITRARY SHIFTS
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

function [] = example_seismic_injection()
nx = 251;
nz = 100;
dx = 10;
xx = (0:nx-1)*dx;
zz = (0:nz-1)*dx;
dt = 0.001;
T  = 0:dt:0.7;
vp  = 1500;

% --- Wavelet
fc = 20;
ricker = @(fm,t) (1-2*pi^2*fm^2*t.^2) .* exp(-pi^2*fm^2.*t.^2);
fs    = ricker(fc, T-0.1 );

% --- Create CENTERED, STAGGERED finite-difference stencil (order = L)
L = 7;
Xur= [-L:L-1];
c = fdweights(0,Xur+1/2,1);
Xu = [L:nx-L]+1; % Updateable weights
Zu = [L:nz-L]+1;

for var=1:7
    % --- Source location
    switch var
        case 1
            sx = 1200; % m
        otherwise
            sx = 1205; % m
    end
    
    % --- Put Receiver 500m further
    rx = sx+500;
    
    % --- Source & receiver stencil
    %     --- Receiver parameters
    N = 7;
    Si = -N:N;
    [ ~, rdx ] = min( abs( xx-rx ) );
    shiftr = (rx-xx(rdx))/dx; % Normalized
    %     --- Source parameters
    M = 4;
    Xi = -M:M;
    [ ~, idx ] = min( abs( xx-sx ) );
    shifts = (sx-xx(idx))/dx; % Normalized
    switch var
        case 1
            a = fdweights(0,Si-shiftr,0);
            b = fdweights(0,Xi-shifts,0);
        case 2
            a = fdweights(0,Si-shiftr,0);
            b = fdweights(0,Xi-shifts,0);
        case 3
            a=FD_rel(Si-shiftr, 0, pi/2,1e4,1e4,false);
            b=FD_rel(Xi-shifts, 0, pi/2,1e4,1e4,false);
        case 4
            a=FD_gv(Si-shiftr, 0, pi/2,1e4,1e4,false);
            b=FD_gv(Xi-shifts, 0, pi/2,1e4,1e4,false);
        case 5
            a=FD_rel(Si-shiftr, 0, pi/2,0,1e4,false);
            b=FD_rel(Xi-shifts, 0, pi/2,0,1e4,false);
        case 6
            a=FD_gv(Si-shiftr, 0, pi/2,0,1e4,false);
            b=FD_gv(Xi-shifts, 0, pi/2,0,1e4,false);
        case 7
            a=FD_hicks( Si-shiftr, N, 10.95 );
            b=FD_hicks( Xi-shifts, M, 6.31 );
    end

    % --- Setup wave system
    p = zeros(nz,nx);
    vx = zeros(nz,nx);
    vz = zeros(nz,nx);

    % --- Recording
    rec = zeros(length(T),1);
    
    % --- Init display
    figure(1)
    clf
    
    % --- Source/receiver arrays
    SRC_c = b'*b;
    REC_c = a'*a;
    
    tic
    for t = 1:length(T)
        % Record P (first sample, like this, is before simulation starts!)
        rec(t) = sum(sum( p(50+Si,rdx+Si).*REC_c )); % A
        
        % Inject p
        p(20+Xi,idx+Xi) = p(20+Xi,idx+Xi) + 2*dt*vp/dx^2*fs(t)*SRC_c * 1e3; % B
        
        % Update P & V
        p(Zu,Xu)   = p(Zu,Xu)    - dt*vp^2*( sum(reshape(c,1,1,length(c)) .* reshape(vx(Zu,Xur+Xu') ,length(Zu),length(Xu),length(c)),3) ) /dx; % C
        p(Zu,Xu)   = p(Zu,Xu)    - dt*vp^2*( sum(reshape(c,1,1,length(c)) .* reshape(vz(Zu'+Xur,Xu)',length(Xu),length(Zu),length(c)),3) )'/dx; % C
        vx(Zu,Xu-1) = vx(Zu,Xu-1)- dt*     ( sum(reshape(c,1,1,length(c)) .* reshape(p( Zu,Xur+Xu') ,length(Zu),length(Xu),length(c)),3) ) /dx; % C
        vz(Zu-1,Xu) = vz(Zu-1,Xu)- dt*     ( sum(reshape(c,1,1,length(c)) .* reshape(p( Zu'+Xur,Xu)',length(Xu),length(Zu),length(c)),3) )'/dx'; % C
        
        imagesc(xx,zz,p)
        axis equal
            hold on
            plot(xx(rdx)+shiftr*dx,zz(50)+shiftr*dx,'o')
            plot(xx(idx)+shifts*dx,zz(20)+shifts*dx,'o')
            hold off
            title(['Case ',num2str(var)','/7'])
            drawnow
    end
    computational_time(var) = toc;

    
    
    % --- Record all simulations
    recs(:,var) = rec;
end

% --- Post-process
[~,inde]=min(abs( T - 0.54) );
RMSE = rms( recs(1:inde,2:end) - recs(1:inde,1) )';% ['Taylor';'rel2';'gv2','relinf','gvinf','hicks'] )
MAXE = max( abs(recs(1:inde,2:end) - recs(1:inde,1)) )

disp(['Taylor rms=',num2str(RMSE(1)),' max error=',num2str(MAXE(1))])
disp(['Hicks  rms=',num2str(RMSE(6)),' max error=',num2str(MAXE(6))])
disp(['Rel l2 rms=',num2str(RMSE(2)),' max error=',num2str(MAXE(2))])
disp(['Relinf rms=',num2str(RMSE(4)),' max error=',num2str(MAXE(4))])
disp(['GV  l2 rms=',num2str(RMSE(3)),' max error=',num2str(MAXE(3))])
disp(['GV inf rms=',num2str(RMSE(5)),' max error=',num2str(MAXE(5))])
fprintf('\n\n')
disp(['The average simulation time = ',num2str(mean(computational_time)),'seconds'])

%%
figure(8)
subplot(2,1,1)
plot( T, recs(:,[1 2 7]), T, recs(:,[3 5 4 6]), '--')
xlim([0.4 0.54])
legend('Reference (on-grid)','Taylor weights','Hicks','|E_{r0}|_2-optimal', ...
    '|E_{r0}|_\infty','|E_{g0}|_2','|E_{g0}|_\infty','Location','NorthWest')
xlabel('Time (s)')
ylabel('Pressure (Pa)')

subplot(2,1,2)

% interpolate finely (sinc interpolation)
recsfft = interpft(recs,50*length(T));
Tnew = linspace( T(1), T(end)+(1-1/50)*dt, length(recsfft) );

plot( Tnew, recsfft(:,[1 2 7])  ,Tnew, recsfft(:,[3 5 4 6]), '--')
xlim([0.4835 0.4848])
ylim([4.78 4.84])
title('Zoomed section')
xlabel('Time (s)')
ylabel('Pressure (Pa)')
%%
end