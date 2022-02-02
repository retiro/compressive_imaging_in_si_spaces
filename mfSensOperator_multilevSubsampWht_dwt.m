function [ y ] = mfSensOperator_multilevSubsampWht_dwt( x, s, mode, picks, n, r_sa, varargin )
% Matrix free sensing operator
%
% Phi - multilevel subsampling, Walsh-Hadamard transform
% Psi - discrete wavelet transform
%
%  Usage
%  -----
%
%   y = mfSensOperator_randSubsampWhtGlobalRand_dwt( x, s, mode, picks, n, r_sa, 'wname' )
%   or
%   y = mfSensOperator_randSubsampWhtGlobalRand_dwt( x, s, mode, picks, n, r_sa, LoR, HiR)
%
%  Input
%  -----
%   x: decomposition vector
%   s: corresponding bookkeeping matrix
%   mode: (1) - forward operator or (2) - transpose
%   picks: random subset of rows of the full measurement matrix Phi
%   n: vector size of the input image (NxN)
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

switch mode
    case 1
        % wavelet reconstruction
        D = waverec2(x,s,LoR,HiR); % B-spline expansion coefficients

        % calculate shift-invariant samples
        SIsamples = conv2(r_sa',r_sa,D,'valid'); % separable 2D convolution
        SIsamples = SIsamples(:); % 2D SIsamples to 1D vector

        % apply fast Walsh-Hadamard transform
        y_all = fastWht(SIsamples')';

        % random subsampling
        y = y_all(picks);
        
    case 2
        % allocate measurements in a full measurement vector
        y_all = zeros(n, 1);
        y_all(picks) = x;
        
        % apply inverse fast Walsh-Hadamard transform
        SIsamples = fastWht(y_all')';
        
        % transposed convolution matrix
        DT = conv2(r_sa',r_sa,reshape(SIsamples, [sqrt(n) sqrt(n)]),'full'); % separable 2D convolution
        
        % wavelet decomposition
        y = wavedec2(DT,size(s,1)-2, wrev(LoR), wrev(HiR))';
        
    otherwise
         error('mode must be 1 (for A*x) or 2 (for A''*x).');
         
end

return