function alpha = calcstepsize(x, d, maxstepsize)
%%PMAT Calculates a stepsize such that the components of the new iteration
%      remain positive.
%
% Copyright (c) 2013.  Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich.

stepsizes = -x./d;
idx_pos = find(stepsizes > 0);
if isempty(idx_pos) 
    alpha = maxstepsize;
else
    alpha = min([0.999*stepsizes(idx_pos); maxstepsize]);
end
          
end
