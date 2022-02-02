function [r_sa] = crossCorr(p_in, p_samp)
%Calculates  the cross-correlation sequence r_sa between B-spline basis functions
%
%  Usage
%  -----
%   r_sa = crossCorr(p_in, p_samp)
%
%  Input
%  -----
%   p_in: order of the B-spline of the generator phi(t)
%   p_samp: order of the B-spline of the sampling kernel zeta(t)
%
%  Output
%  ------
%   r_sa: calculated cross-correlation sequence
%
% Copyright (c) 2021,  Tin Vlasic
% ---------------------------------------------------------------
syms t;
% generate signal basis functions
N_overlap_sig = floor((p_in + p_samp + 1)/2); % number of overlaping basis functions
beta = 0;
for j = 0 : p_in + 1
   beta = beta + ((-1)^j)/(factorial(p_in))*nchoosek(p_in+1,j)* ... 
                  (t-1/2+(p_in+1)/2-j)^p_in*heaviside(t-1/2+(p_in+1)/2-j);
end
B_a = sym(zeros(1, 1+2*N_overlap_sig));
for i = -N_overlap_sig : N_overlap_sig
    B_a(i+N_overlap_sig+1) = subs(beta, t, t-i);
end

% generate sampling basis functions
N_overlap_samp = floor((p_samp+1)/2); % number of overlaping basis functions
beta = 0;
for j = 0 : p_samp+1
   beta = beta + ((-1)^j)/(factorial(p_samp))*nchoosek(p_samp+1,j)* ... 
                  (t-1/2+(p_samp+1)/2-j)^p_samp*heaviside(t-1/2+(p_samp+1)/2-j);
end
B_s = sym(zeros(1, 1+2*N_overlap_samp));
for i = -N_overlap_samp : N_overlap_samp
    B_s(i+N_overlap_samp+1) = subs(beta, t, t-i);
end

% calculate cross-correlation sequence
r_sa = zeros(1, 2*N_overlap_sig+1);
for i = 1 : 2*N_overlap_sig + 1
    r_sa(i)= int(B_s(N_overlap_samp + 1)*B_a(i), -(N_overlap_samp+1), (N_overlap_samp+1));
end

% cross-correlation sequence
cross_vec = sprintf('%1.4f ', r_sa);
fprintf('Cross-correlation sequence: r_sa=[ %s] \n', cross_vec);

end

