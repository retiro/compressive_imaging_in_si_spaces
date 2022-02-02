function x = Pmat(T, m, n, y)
%PMAT Implements the preconditioner 
%
%  P^{-1} = [I, rho*(rho*I+T_1)^{-1}]*[(rho*I+T_1)^{-1}, 0]*[         I          , 0]
%           [0,          I          ] [       0        , S] [(rho*I+T_1)^{-1}*rho, I]
%
%  where S = rho*I+T_2 - rho^2*(rho*I+T_1)^{-1}.
%
%  for the matrix
%
%  H = [ A^T*A, -A^T*A] + [T_1,  0 ]
%      [-A^T*A,  A^T*A]   [ 0 , T_2]
%
% Copyright (c) 2012.  Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich.

x = y;
rho = m/n;
D1 = T(1:n) + rho;
D2 = T(n+1:end) + rho;

x(n+1:2*n) = x(n+1:2*n) + rho*x(1:n)./D1;

x(1:n) = x(1:n)./D1;
x(n+1:2*n) = x(n+1:2*n)./(D2-rho^2./D1);

x(1:n) = x(1:n) + rho*x(n+1:2*n)./D1;

end
