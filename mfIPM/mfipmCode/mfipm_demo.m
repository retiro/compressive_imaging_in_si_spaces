%MFIPM_DEMO Demonstrates the usage of the Matrix-free IPM solver.
%  
%  Two simple examples are shown.
%    1) One dimensional real Compressed Sensing problem with implicitily
%       stored Gaussian sensing matrix.
%    2) One dimensional real Compressed Sensing problem using an 1-d
%       partial DCT operator as a sensing matrix in the form of an A_operator.
%       ** Requires Signal Processing Toolbox **
%
% Copyright (c) 2013.  Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich.

% ------------------------------------------------------
% Example 1, implicitily stored Gaussian sensing matrix.
% ------------------------------------------------------
% fprintf('Example 1, Gaussian sensing matrix \n \n');

% Fix the seed.
randn('state',200);
rand('state',200);

% Number of measurements.
m = 2^17;

% Size of the signal.
n = 2^18;

% Number of nonzeros in the sparsest representation.
k = 2^15;

% % Create the sensing matrix. Here we assume that the sparsity matrix
% % is the identity.
% A = randn(m,n);
% 
% % Create artificially an exact sparse representation.
% p     = randperm(n);
% x_opt = zeros(n,1); 
% 
% x_opt(p(1:k)) = 2*randn(k,1);
% 
% % Calculate the measurements.
% b = A*x_opt;
% 
% % Set the regularization parameter.
% tau = 1.0e-3;
% 
% % Options.
% opts.tol       = 1.0e-10;
% opts.maxiters  = 100;
% opts.tolpcg    = 1.0e-2;
% opts.mxiterpcg = 200;
% opts.verbose   = 2;
% 
% % Call the Matrix-free IPM solver.
% [x out] = mfipm(n, A, b, tau, opts);
% 
% % Print info.
% fprintf(['\n||x-x_opt||/||x_opt|| = %1.3e, Iters = %3d,'               ,...
%          ' MatVecs = %5d \n \n'], norm(x-x_opt)/norm(x_opt), out.iters ,...
%          out.nMat);

% ---------------------------------------------------------------
% Example 2, 1-d DCT sensing matrix in the form of an A_operator.
% ---------------------------------------------------------------
fprintf('Example 2, partial DCT matrix, ** Requires Signal Processing Toolbox ** \n \n');

% Partial 1-d discrete cosine transform. Here we assume that the sparsity
% matrix is the identity.
picks = randperm(n);
A = A_operator(@(x) pdct(x,1,n,sort(picks(1:m))),...
               @(x) pdct(x,2,n,sort(picks(1:m)))); 

% Create artificially an exact sparse representation.
p     = randperm(n);
x_opt = zeros(n,1); 

x_opt(p(1:k)) = 2*randn(k,1);

% Calculate the measurements.
b = A*x_opt;

tau = 1.0e-4; % regularization parameter.

% Options.
opts.tol       = 1.0e-10;
opts.maxiters  = 100;
opts.tolpcg    = 1.0e-4;
opts.mxiterpcg = 200;
opts.verbose   = 1;

% Call the Matrix-free IPM solver.
[x out] = mfipm(n, A, b, tau, opts);

% Print info.
fprintf(['\n||x-x_opt||/||x_opt|| = %1.3e, Iters = %3d,'               ,...
         ' MatVecs = %5d \n \n'], norm(x-x_opt)/norm(x_opt), out.iters ,...
         out.nMat);

figure; stem(x(1:100))
hold on; stem(x_opt(1:100))

%clear;
