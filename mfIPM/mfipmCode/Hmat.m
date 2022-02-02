function y = Hmat(A, x, T)
%HMAT Performs matrix-vector products with the following matrix
%
%  [ A^T*A, -A^T*A] + [T_1,  0 ]
%  [-A^T*A,  A^T*A]   [ 0 , T_2]
%
% Copyright (c) 2012.  Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich

n2 = length(T);
n  = n2/2;

y = zeros(n2,1);

temp = A'*(A*(x(1:n)-x(n+1:end)));

y(1:n) = T(1:n).*x(1:n) + temp;
y(n+1:end) = T(n+1:end).*x(n+1:end) - temp;

end
