function [ y ] = mfSensOperator_randSubsampWhtGlobalRand_dct( x, mode, picks, perm, permReturn, r_sa )
% Matrix free sensing operator
%
% Phi - random subsampling, Walsh-Hadamard transform, global randomizer
% Psi - discrete cosine transform
%
%  Usage
%  -----
%
%   y = mfSensOperator_randSubsampWhtGlobalRand_dct( x, mode, picks, perm, permReturn, r_sa)
%
%  Input
%  -----
%   x: decomposition vector
%   mode: (1) - forward operator or (2) - transpose
%   picks: random subset of rows of the full measurement matrix Phi
%   perm: global randomizer locations (scrambles the signal's sample locations globally)
%   permReturn: indices for reverse permutation
%   r_sa: cross-correlation sequence between sampling kernel and generator in 1D
%                        
%  Output
%  ------
%   y: mode (1) - y = A*x
%      mode (2) - y = A'*x
%
% Copyright (c) 2021,  Tin Vlasic and Damir Sersic
% ---------------------------------------------------------------
n = numel(perm); % coefficients length
x_len_1d = sqrt(numel(x));

switch mode
    case 1
        % wavelet reconstruction
        x = reshape(x, [x_len_1d, x_len_1d]);
        D = idct2(x); % B-spline expansion coefficients
        
        % calculate shift-invariant samples
        SIsamples = conv2(r_sa',r_sa,D,'valid'); % separable 2D convolution
        SIsamples = SIsamples(:); % 2D SIsamples to 1D vector

        % measurement with structurally random matrix
        % globally randomize SI sample locations
        SIsamples_perm = SIsamples(perm);

        % apply fast Walsh-Hadamard transform
        y_all = fastWht(SIsamples_perm')';

        % random subsampling
        y = y_all(picks);
        
    case 2
        % allocate measurements in a full measurement vector
        y_all = zeros(n, 1);
        y_all(picks) = x;
        
        % apply inverse fast Walsh-Hadamard transform
        SIsamples_perm = fastWht(y_all')';
        
        % return measurements to its original positions
        SIsamples = SIsamples_perm(permReturn);
        
        % transposed convolution matrix
        DT = conv2(r_sa',r_sa,reshape(SIsamples, [sqrt(n) sqrt(n)]),'full'); % separable 2D convolution
        
        % wavelet decomposition
        y = dct2(DT);
        y = y(:);
        
    otherwise
         error('mode must be 1 (for A*x) or 2 (for A''*x).');
         
end

return