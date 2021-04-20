function [diff]=eldiff(x)

% Computes an array of differences between successive elements of input matrix
%
% INPUTS: 
%   x: 1D of length n
%
% OUTPUTS:
%   diff: 1D array of length n-1 where each element represents the
%         difference between x(i) and x(i+1)
%
% Written by Danica Roth, University of Oregon, January 2017.
%

n=length(x);

diff=zeros((n-1),1);
for i=1:(n-1)
    diff(i)=x(i+1)-x(i);
end