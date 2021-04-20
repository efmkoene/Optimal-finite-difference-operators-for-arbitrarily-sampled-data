%% Script to recreate the Hicks (2002) coefficients
% from 'Arbitrary source and receiver positioning in finite-difference
% schemes using Kaiser windowed sinc functions'

% from 'Optimal finite-difference operators for arbitrarily sampled data'
% Copyright 2020, SEG. Erik Koene & Johan Robertsson



function [ weights ] = FD_hicks( x,r,b )
%FD_HICKS computes the weights for Hicks' proposed monopole injection scheme
% Input x as unit distance properties
% Example:
%   FD_hicks( [-5:4]+0.5,4,6.81 )

kaiserwin =@(x,r,b) besseli(0,b*sqrt(1-(x/r).^2)) ./ besseli(0,b) .* (abs(x)<r);
weights = sinc( x ) .* kaiserwin( x, r, b );

end

