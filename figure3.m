%% Script to recreate Figure 3,
% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson

clear all;

% --- Reset random number generator
rng('default');

% --- Medium extent
x0  = 0;
xend= 3000;
N = 301;

% --- Create the grid (eq. 34)
ss_cos = pi/6;
x1 = cos(ss_cos) - cos( ss_cos + [0:N-1] / (N-1) * (pi-2*ss_cos) );
x1 = x1 / (2*cos(ss_cos));
x1 = x0 + x1*xend;
x = x1 + (randn(1,N));
x = sort(x);
x(end) = 3000-x(1);

% --- Create the figure
figure(3)
subplot(5,1,1)
plot( x, x*0, 'k.' )
title('Node positions x (left edge)')
xlabel('x (m)')
xlim([0 500])

subplot(6,1,[3:6])
L = 6
[~,idx]=min( abs( x - 500 ) );
for k=1:idx-1
    % All the internal points
    if k>L && (N-k)>L
        smp = k-[-L:L];
    elseif k<=L
        smp = [1:2*L+1];
    elseif (N-k)<=L
        smp = [N-2*L:N];
    end
    plot( x(1:idx) , k, '.', 'Color',[0 0 0]+0.8)
    hold on
    plot( x(smp), k, 'r.' )
    plot( x(k), k, 'k.')
end
xlabel('x (m)')
ylabel('Node i')
title('Stencils for node positions x_i (left edge)')
xlim([-1 500])
ylim([0 56])

