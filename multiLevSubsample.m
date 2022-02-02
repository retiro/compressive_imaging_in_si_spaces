function [ mLevPicks ] = multiLevSubsample( n, I, M )
%multiLevSubsample Subsamples measurement matrix in levels
% For asymptotically incoherent CS matrices and asymptotically sparse signals
%
%  Usage
%  -----
%   mLevPicks = multiLevSubsample( n, I, M )
%
%  Input
%  -----
%   n: measurement mask resolution (NxN)
%   I: number of scale levels
%   M: number of measurements
%
%  Output
%  ------
%   mLevPicks: multilevel subset of rows of the full measurement matrix Phi
%
% Copyright (c) 2021,  Tin Vlasic
% ---------------------------------------------------------------
%%
lev = I+2;  % number of subsampling levels
mLevPicks = cell(1, lev); % multilevel picks

% coefficient percentage in a decomposition
p = zeros(1, lev);
p(1) = 1/2^2;
for i = 2:I
    p(i) = 2 ./ (2^(2*(i-1))) + 1/(2^2^i);
end
p(lev-1) = 2/(2^2^I);
p(lev) = 1/(2^2^I);

% number of samples in each level
w = (1 : lev)/lev; % weights
a = M / sum(w .* p);
m = a .* w .* p; % number of picks in each level
mMax = p*n; % max number of picks in each level

for i = lev:-1:1
    % calculate residual and transfer it to the next level
    res = round(m(i)) - mMax(i);
    if res > 0 
        try m(i-1) = m(i-1) + res; catch; end
        m(i) = mMax(i);
    else
        m(i) = round(m(i));
    end
    % random subsampling in levels
    mLevPicks{1,i} = sort(randsample( 1+sum(mMax(lev:-1:i+1)) : sum(mMax(lev:-1:i)), m(i)));
end

end

