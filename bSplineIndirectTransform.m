function [sig] = bSplineIndirectTransform(coeffs, n, m)
%Reconstructs signal from B-spline coefficients
%
%  Usage
%  -----
%   sig = bSplineIndirectTransform(coeffs, n, m)
%
%  Input
%  -----
%   coeffs: B-spline coefficients
%   n: order of the B-spline kernel
%   m: upsampling factor
%
%  Output
%  ------
%   sig: pointwise values of the signal
%
% Copyright (c) 2021,  Tin Vlasic
% ---------------------------------------------------------------
%% discrete B-spline basis function
b0_1 = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]; % starting values
c0_1 = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]; % starting values

lng = length(b0_1);
bn_1 = b0_1; % allocate bn_1
cn_1 = c0_1; % allocate cn_1
for i = 1 : n
    for k = 1+1 : lng-1
        bn_1(k) = ( ((k-(lng+1)/2)+(i+1)/2)*c0_1(k) + ...
                    (-(k-(lng+1)/2)+(i+1)/2)*c0_1(k-1)) / (i);
    end
    for k = 1+1 : lng-1
        cn_1(k) = ( ((k-(lng+1)/2)+(i+2)/2)*b0_1(k+1) + ...
                     (-(k-(lng+1)/2)+i/2)*b0_1(k)  ) / (i);
    end
    b0_1 = bn_1;
    c0_1 = cn_1;
end

%% discrete B-spline basis function with upsampling integer m
b0_m = [ ones(1, floor(m/2)), 1, ones(1, floor((m-1)/2)) ]; % starting values

bn_m = b0_m;
for i = 1 : n
    bn_m = conv(bn_m, b0_m);
end

if mod(m,2) ~= 0
    % m odd
    bn_m = 1/(m^n) * conv(bn_m, bn_1);
elseif mod(n,2) ~= 0
    % n odd and m even
    bn_m = 1/(m^n) * conv(bn_m, bn_1);
else
    % n even and m even
    bn_m = 1/(m^n) * conv(bn_m, cn_1);
end
%% calculate indirect B-spline filter's transfer function
num = nonzeros(bn_m)'; % take only nonzero values from b^n_m

%% upsamle the coefficients if needed
d_up = zeros([size(coeffs,1)*m 1]);
d_up(1:m:end,1:m:end) = coeffs;

%% filter calculated B-spline coeffiecients with the indirect B-spline filter
if n ~= 0
    N_delay = (length(num)-1)/2; % delay
else 
    N_delay = 0;
end
N_overlap = floor((n+1)/2);
n_tilde = floor(n/2);

sig = filter(num, 1, d_up); % causal filtering
if mod(n,2) == 0
    % n even
    sig = sig(1 + m*n : end);
else
    % n odd
    sig = sig(1 + N_delay + m * N_overlap : end - (N_overlap - n_tilde) );
end
