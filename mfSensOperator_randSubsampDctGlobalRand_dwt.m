function [ y ] = mfSensOperator_randSubsampDctGlobalRand_dwt( x, s, mode, picks, perm, permReturn, r_sa, varargin )
% Matrix free sensing operator
%
% Phi - random subsampling, discrete cosine transform, global randomizer
% Psi - discrete wavelet transform
%
%  Usage
%  -----
%
%   y = mfSensOperator_randSubsampDctGlobalRand_dwt( x, s, mode, picks, perm, permReturn, r_sa, 'wname')
%   or
%   y = mfSensOperator_randSubsampDctGlobalRand_dwt( x, s, mode, picks, perm, permReturn, r_sa, LoR, HiR)
%
%  Input
%  -----
%   x: decomposition vector
%   s: corresponding bookkeeping matrix
%   mode: (1) - forward operator or (2) - transpose
%   picks: random subset of rows of the full measurement matrix Phi
%   perm: global randomizer locations (scrambles the signal's sample locations globally)
%   permReturn: indices for reverse permutation
%   r_sa: cross-correlation sequence between sampling kernel and generator in 1D
%   varargin: 'wname' - a character vector containing the wavelet name
%             LoR, HiR - instead of giving the wavelet name, you can give the filters
%                        
%  Output
%  ------
%   y: mode (1) - y = A*x
%      mode (2) - y = A'*x
%
% Copyright (c) 2021,  Tin Vlasic and Damir Sersic
% ---------------------------------------------------------------
%narginchk(7,8)
if ischar(varargin{1})
    [LoR, HiR] = wfilters(varargin{1}, 'r');
else
    LoR = varargin{1}; HiR = varargin{2};
end
n = numel(perm); % coefficients length

switch mode
    case 1
        % wavelet reconstruction
        D = waverec2(x,s,LoR,HiR); % B-spline expansion coefficients

        % calculate shift-invariant samples
        SIsamples = conv2(r_sa',r_sa,D,'valid'); % separable 2D convolution
        SIsamples = SIsamples(:); % 2D SIsamples to 1D vector

        % measurement with structurally random matrix
        % globally randomize SI sample locations
        SIsamples_perm = SIsamples(perm);

        % apply fast Walsh-Hadamard transform
        y_all = dct(SIsamples_perm);

        % random subsampling
        y = y_all(picks);
        
    case 2
        % allocate measurements in a full measurement vector
        y_all = zeros(n, 1);
        y_all(picks) = x;
        
        % apply inverse fast Walsh-Hadamard transform
        SIsamples_perm = idct(y_all);
        
        % return measurements to its original positions
        SIsamples = SIsamples_perm(permReturn);
        
        % transposed convolution matrix
        DT = conv2(r_sa',r_sa,reshape(SIsamples, [sqrt(n) sqrt(n)]),'full'); % separable 2D convolution
        
        % wavelet decomposition
        y = wavedec2(DT,size(s,1)-2, wrev(LoR), wrev(HiR))';
        
    otherwise
         error('mode must be 1 (for A*x) or 2 (for A''*x).');
         
end

return