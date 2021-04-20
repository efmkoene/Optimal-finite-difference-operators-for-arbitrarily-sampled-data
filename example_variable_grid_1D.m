%% Do FD modeling with variable grid spacing, 1D. I will create a continuously varying grid size.
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

function [] = example_variable_grid_1D()
%%
x0  = 0;
xend= 3000;
N = 301;

% Create the variable grid space (with fixed random number generator)
% --- Get subset of the cosine (from ss_cos : pi-ss_cos, rather than from 0 to pi)
ss_cos = pi/6;
x1 = - cos( ss_cos + [0:N-1] / (N-1) * (pi-2*ss_cos) ) + cos(ss_cos);
x1 = x1 / (2*cos(ss_cos));
x1 = x0 + x1*xend;
rng('default');
x = x1 + (randn(1,N))*1;
x = sort(x);
x(end) = 3000-x(1);

% --- Init
sx = 75; % Inject source in ..th index
rx = 150;
c = 2000 * ones(N,1);
% --- Time parameters
dt = 0.001;
T  = 0:dt:5;
nt = length(T);
record = zeros(nt,3);

% --- pre-Compute Forward Time Dispersion Transform, cf. Koene et al. (2018); Eliminating time dispersion from seismic wave modeling; GJI
A= exp( -1i * 2*sin([0:nt]*pi/(2*nt))' * [0:2*nt-1] );
B = ifft( A, 2*nt, 'symmetric');
FTDT = B(1:nt,1:nt);

% --- pre-Compute Inverse Time Dispersion Transform, cf. Koene et al. (2018); Eliminating time dispersion from seismic wave modeling; GJI
A1 = exp( -1i * 2 * sin([0:nt-1]*pi/(2*nt))' * [0:2*nt-1] );
v  = cos([0:nt-1]*pi/(2*nt));
A2 = v' .* A1;
B  = ifft( A2, 2*nt, 'symmetric' );
ITDT  = B(1:nt,1:nt)';

% --- Wavelet
fc = 30;          % Choose 25 Hz or 50 Hz (or anything else!)
ricker = @(fm,t) (1-2*pi^2*fm^2*t.^2) .* exp(-pi^2*fm^2.*t.^2);
fs0= ricker(fc, T-1.5/fc );
fs = FTDT * fs0(:);   

% --- Wavefield recording
record_everything=zeros(N,nt);


for L = [6 14 22]
    
    % --- Choose the range for which the median 'dx' is established
    dxf=L;

% Create N stencils for the variable grid
stn = zeros(N,2*L+1);
tic
for k=1:N
    % All the internal points
    if k>L && (N-k)>L
        dxsmp = k-[-dxf:dxf];
        smp = k-[-L:L];
%         stn(k,:) = fdweights( 0, x(k)-x(smp), 2 );
        fact = 1.25;
    elseif k<=L
        dxsmp = [1:2*dxf+1];
        smp = [1:2*L+1];
%         stn(k,:) = fdweights( 0, x(k)-x(smp), 2 );
        fact = 1;
    elseif (N-k)<=L
        dxsmp = [N-2*dxf:N];
        smp = [N-2*L:N];
%         stn(k,:) = fdweights( 0, x(k)-x(smp), 2 );
        fact = 1;
    end
    
    % --- Median dx
    dx = median(abs(diff( x(k)-x(dxsmp) ))) * fact;
    % --- Least-squares relative optimal coefficients. Suppress visual output
    stn(k,:) = FD_rel( x(k)-x(smp), 2, pi/dx, 1e3, [], false );

end
stencil_time(floor(L/6)) = toc;



% Space parameters
p_new = zeros(N,1);
p     = zeros(N,1);
p_old = zeros(N,1);

figure(1)
clf
tic
for t=1:nt
    
    % Record data
    record(t,floor(L/6)) = p(rx);
   
    % Source injection
    p(sx) = p(sx) + fs(t)*dt^2*c(sx)^2;
    
    % Update equation
        for k=1:N
            if k>L && (N-k)>L
                samp = p(k-[-L:L]');
            elseif k<=L
                samp = p([1:2*L+1]');
            elseif (N-k)<=L
                samp = p([N-2*L:N]');
            end
            % The update equation
            p_new(k) = c(k)^2*dt^2 * (stn(k,:) * samp) + 2*p(k) - p_old(k);
        end
    
     % Boundary conditions (free surface) >> Do not update smallest elements
     p_new(1) = 0;
     p_new(N) = 0;
     
     % Prepare for next update
     p_old = p;
     p = p_new;
     
     % Record full wavefield
     if L==22
         record_everything(:,t) = p(:);
     end
    
     if mod(t,5)==0
     plot( x, p , '.-' )
     title(['t=',num2str(dt*t)])
     xlabel('x (m)')
     drawnow
     end
    
end
computational_time(floor(L/6)) = toc;

end

for L=1:3
    rect(:,L) = ITDT*record(:,L);
end


% --- Create the output figure
figure(4)
subplot(2,1,1)

% --- Filtere the whole dataset
record_ITDT = record_everything;
for i=1:N
    record_ITDT(i,:) = ITDT * record_everything(i,:)';
end

pcolor( x, T, record_ITDT' )
caxis(0.1*caxis)
ylabel('Time (s)')
xlabel('Offset x (m)')
shading interp
axis ij tight
title('1D wave propagation in space-time domain')

subplot(2,1,2)

len = 1:1400;

plot( ... % The L=6 case
      T(len), rect(len,1), 'b', ...
      T(len), rect(len+3000,1)+85, 'k', ...
      T(len), 50*(rect(len,1) - rect(len+2999,1))+170, 'r', ...
      ... % The L=14 case
      T(len), rect(len,2)+350, 'b', ...
      T(len), rect(len+3000,2)+435, 'k', ...
      T(len), 50*(rect(len,2) - rect(len+2999,2))+520, 'r', ...
      ... % The L=22 case
      T(len), rect(len,3)+650, 'b', ...
      T(len), rect(len+3000,3)+735, 'k', ...
      T(len), 50*(rect(len,3) - rect(len+2999,3))+820, 'r')
ylim([-150 900])
xlim([0 1.4])
xlabel('t_0 + t (s)')
view([90 90])
text(0,85,'L=6')
text(0,435,'L=14')
text(0,735,'L=22')

disp(['Stencil-obtaining time=',num2str(stencil_time(1)),'. Simulation time = ',num2str(computational_time(1)),'seconds for L=6'])
disp(['Stencil-obtaining time=',num2str(stencil_time(2)),'. Simulation time = ',num2str(computational_time(2)),'seconds for L=14'])
disp(['Stencil-obtaining time=',num2str(stencil_time(3)),'. Simulation time = ',num2str(computational_time(3)),'seconds for L=22'])

end