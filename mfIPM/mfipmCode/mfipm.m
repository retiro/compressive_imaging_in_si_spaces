function [x out] = mfipm(n, A, b, tau, opts)
%MFIPM Matrix-free Interior Point solver.
%
%  The Matrix-free IPM solver for one dimensional Compressed Sensing
%  problems proposed in [1].
%
%    [1] Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich 
%        "Matrix-Free Interior Point Method for CS problems", 
%        ERGO Technical Report, 2012.
%
%  The Matrix-free IPM solves the l1-regularised program:
%
%    min tau*||x||_1 + 0.5*||Ax - b||^2,
%
%  using an infeasible primal-dual interior point algorithm, for details
%  please see at [1].
%
%  Usage
%  -----
%
%    [x out] = mfipm(n, A, b, tau, opts) 
%
%  Input (required)
%  ----------------
%
%    n:   Number of columns of the constraint matrix A.
%    A:   The constraint matrix.
%    b:   The right hand side of the constraints Ax=b.
%    tau: The regularisation parameter.
%
%  Input (optional)
%  ----------------
%  
%    The opts argument is a structure and it is optional. Here is a 
%    short explanation of the available options for opts:
%
%    opts.tol:       Optimality tolerance of the Matrix-free IPM.
%                    Default value: 1.0e-8.
%    opts.maxiters:  Maximum iterations of the Matrix-free IPM.
%                    Default value: 100.
%    opts.tolpcg:    Optimality tolerance for the PCG method.
%                    Default value: 1.0e-4.
%    opts.mxiterpcg: Maximum iterations of the PCG method.
%                    Default value: 200.
%    opts.maxstepsz: Maximum stepsize for primal and dual directions.
%                    Default value 0.999.
%    opts.cnt1:      Centering parameter when predictor directions are 
%                    performed. Set to 0 for pure predictor directions.
%                    Default value 0.1.
%    opts.cnt2:      Centering parameter when predictor directions are 
%                    performed and one of the primal and dual stepsizes of 
%                    the previous iteration were less than opts.ap for last
%                    iteration. 
%                    Set to 0 for pure predictor directions.
%                    Default value 0.5.
%    opts.cnt3:      Centering parameter when corrector directions are 
%                    performed. Set to 1 for perfect centering to the
%                    central path.
%                    Default value 0.8.
%    opts.ap:        If primal or dual stepsize <= ap then use cnt2 as as
%                    centering parameter, otherwise use cnt1.
%                    Default value 0.5.
%    opts.ac:        If primal or dual stepsize after calculating the
%                    predictor direction are <= ac then perform a corrector
%                    update.
%                    Default value 0.1.
%    opts.linsol:    0 to use the Conjugate Gradients method as a linear 
%                    system solver, 1 to use a direct solver.
%                    Default value: 0.
%    opts.verbose:   0 for no printing info, 1 for basic printing info
%                    and 2 for detailed printing info.
%                    Default value: 1.
%
%    ***************************************************************************************
%    If the measurements are noiseless please use opts.tol >= 1.e-8. In case that the
%    measurements are contaminated with noise then please use opts.tol >= 1.0e-6. Regarding
%    the accuracy of the PCG method opts.tolpcg=1.0e-1 to 1.0e-4 should be enough.
%    ***************************************************************************************
%
%  Output
%  ------
%
%    x:               The reconstructed sparse representation.
%    out.totcputime:  Total CPU time required by the Matrix-free IPM.
%    out.iters:       Number of Matrix-free IPM iterations.
%    out.cpusteps:    CPU time used for the calculation of the stepsizes.
%    out.totpcgiters: Total number of PCG iterations.
%    out.nMat:        Number of matrix-vector products.
%    out.predictors:  Number of predictor directions.
%    out.correctors:  Number of corrector directions.
%    out.linsystems:  Total number of linear systems solved. It is equal
%                     to out.predictors + out.correctors.
%    out.status:      1 if the problem was solved successfully, 0 if 
%                     maximum number of iterations reached.
%    out.message:     Termination status of the Matrix-free IPM
%                     solver. One of the following two messages is returned:
%                     i)  'Converged'
%                     ii) 'Maximum iterations reached'
%
%
%  Note
%  ----
%
%    If A is an operator, then it should be passed to the Matrix-free IPM 
%    as an A_operator type. To use A_operator variables please add to your 
%    MATLAB (if it is not already added) the subdirectory:
%    
%      "main mfipm directory"/operators
% 
%    Example: To create an A_operator object, two functions or function 
%    handles for computing A*x and A'*x, respectively, must be given. 
%    Suppose they are AA and AT,
%          AA: a function handle such that AA(x) = A*x,
%          AT: a function handle such that AT(x) = A'*x.
%    Then A can be created as an A_operator by
%          A = A_operator(@(x) AA(x), @(x) AT(x));
%
%    An example for A being an implicit matrix that performs a discrete 
%    cosine transform (DCT) and returns the subset of the DCT result 
%    corresponding to omega is provided below
%
%     function y=pdct(x,picks); y=dct(x); y=y(picks); end
%     function y=pidct(x,n,picks); y=zeros(n,1); y(picks)=x; y=idct(y); end
%     A = A_operator(@(x) pdct(x,omega), @(x) pidct(x,n,omega)); 
%
%    **This example has been copied from the FPC_AS.m in the 
%      FPC_AS_v1.21 software package.
%      url: http://www.caam.rice.edu/~optimization/L1/FPC_AS/ 
%
% Copyright (c) 2012.  Kimon Fountoulakis, Jacek Gondzio and Pavel Zhlobich.

% Check number of arguments.
if nargin < 4
    error('At least four arguments are required.');
end

% Check mandatory arguments.
if ~isa(A, 'A_operator') && ~(sum(size(A) > 1) == 2)
    error('The argument "A" must be a matrix or of type "A_operator".');
end
sizes_b = size(b);
if sizes_b(1) > 1 && sizes_b(2) > 2
    error('The argument "b" should be a vector.');
end
if sizes_b(2) > 1
    b = b';
end
if length(b) > n
    error('"m" should be less than "n".');
end

% Check optinal arguments.
opts = chkOptargs(opts);

% Counters.
out.predictors  = 0;
out.correctors  = 0;
out.totpcgiters = 0;
out.linsystems  = 0;
out.nMat        = 0;
out.iters       = 0;
out.cpusteps    = 0;
out.totcputime  = 0;

% Initialize sizes.
n2 = 2*n;
m  = length(b);

% CPU time.
tstart = cputime;

% Calculate linear part coefficients.
ATb      = A'*b;
c        = tau + [-ATb;ATb];
out.nMat = out.nMat + 1;

% Calculate initial point here.
reg  = 1.0e-2*ones(n2,1);
H    = @(q) Hmat(A, q, reg);
Prec = @(p) Pmat(reg, m, n, p);
    
[z,~,~,itersPcg] = pcg(H,-c,1.0e-2,200,Prec);
out.nMat = out.nMat + 2*itersPcg;
    
z = abs(z) + norm(b);
s = norm(b)*ones(n2,1);

% Initialize the rest of the variables here.
ATAblkmat = A'*(A*(z(1:n)-z(n+1:end)));
out.nMat  = out.nMat + 2;

mu        = z'*s/n2;
dual_e    = s - c - [ATAblkmat;-ATAblkmat];
normc     = norm(c);

if opts.linsol
   F   = A'*A;
   FFt = [F,-F;-F,F];
end

% Print info.
if opts.verbose == 2
   fprintf(' %5s | %5s |  %5s   | %8s | %8s | %9s | %10s |%10s \n', ' Info ', 'Iters',...
        'rDGAP', 'rDUALinf','Iters PCG', 'RelRes PCG', 'Primal Step', 'Dual Step');
   fprintf('  -------------------------------------------------------------------------------------- \n');
end

while 1

    % Theta^{-1} diagonal matrix
    T = s./z;
        
    % Compute affine scaling direction.
    if out.iters ~= 0 && (alphaP < 0.5 || alphaD < 0.5)
        sigma_a = opts.cnt2;
    else
        sigma_a = opts.cnt1;
    end
    
    barrier_a = sigma_a*mu;
    
    f2 = barrier_a - z.*s;
    
    invZf2 = f2./z;
    if ~opts.linsol
        H    = @(q) Hmat(A, q, T);
        Prec = @(p) Pmat(T, m, n, p);
    
        [dza,~,rrPcg,itersPcg] = pcg(H,dual_e+invZf2,opts.tolpcg,opts.mxiterpcg,Prec);
        out.nMat = out.nMat + 2*itersPcg;
    else
        H = FFt + diag(T);
        dza = H\(dual_e+invZf2);
    end
    dsa = invZf2 - T.*dza;
    
    % Choose stepsize for the affine scaling direction.
    tstepsstart = cputime; 
    
    alphaP = calcstepsize(z,dza,opts.maxstepsz);
    alphaD = calcstepsize(s,dsa,opts.maxstepsz);
    
    tstepsend = cputime;
    out.cpusteps = out.cpusteps + (tstepsend - tstepsstart);
    
    z_new = z + alphaP*dza;
    s_new = s + alphaD*dsa;
    
    % Print info.
    if opts.verbose == 2 && ~opts.linsol
        fprintf('  %5s | %3d   | %8s | %8s | %5d     |  %1.2e  |   %1.2e  |  %1.2e \n',...
        'Pred.', out.iters, '', '', itersPcg, rrPcg, alphaP, alphaD);
    elseif opts.verbose == 2 && opts.linsol
        fprintf('  %5s | %3d   | %8s | %8s |   %1.2e  |  %1.2e \n',...
        'Pred.', out.iters, '', '', alphaP, alphaD);      
    end
    
    out.predictors  = out.predictors + 1;
    out.totpcgiters = out.totpcgiters + itersPcg;
    
    % Compute centering direction, if needed.
    if alphaP < 0.1 || alphaD < 0.1
        
       % Compute centering parameter.
       ATAblkmat   = A'*(A*(z_new(1:n)-z_new(n+1:end)));
       out.nMat    = out.nMat + 2;
       dual_e_c  = s_new - c - [ATAblkmat;-ATAblkmat];
       sigma_c   = opts.cnt3;
       mu_c      = z_new'*s_new/(n2);
       barrier_c = sigma_c*mu_c;
       
       f2_c = barrier_c - z_new.*s_new;
    
       invZf2_c = f2_c./z_new;
       if ~opts.linsol 
           [dzc,~,rrPcg,itersPcg] = pcg(H,dual_e_c+invZf2_c,opts.tolpcg,...
               opts.mxiterpcg,Prec);
           out.nMat = out.nMat + 2*itersPcg;
       else
           dzc = H\(dual_e_c+invZf2_c);
       end
       dsc = invZf2_c - T.*dzc;
       
       % Compute total direction.
       dz = dza + dzc;
       ds = dsa + dsc;
    
       % Choose stepsize for the final direction.
       tstepsstart = cputime;
       
       alphaP = calcstepsize(z,dz,opts.maxstepsz);
       alphaD = calcstepsize(s,ds,opts.maxstepsz);
       
       tstepsend    = cputime;
       out.cpusteps = out.cpusteps + (tstepsend - tstepsstart);
    
       z_new = z + alphaP*dz;
       s_new = s + alphaD*ds;
       
       % Print info.
       if opts.verbose == 2 && ~opts.linsol
           fprintf('  %5s | %3d   | %8s | %8s | %5d     |  %1.2e  |   %1.2e  |  %1.2e \n',...
           'Corr.', out.iters, '', '', itersPcg, rrPcg, alphaP, alphaD);
       elseif opts.verbose == 2 && opts.linsol
           fprintf('  %5s | %3d   | %8s | %8s |   %1.2e  |  %1.2e \n',...
           'Corr.', out.iters, '', '', itersPcg, rrPcg, alphaP, alphaD);               
       end
       
       out.correctors  = out.correctors + 1;
       out.totpcgiters = out.totpcgiters + itersPcg;
    end
    
    z = z_new;
    s = s_new;

    % Termination criteria.
    ATAblkmat = A'*(A*(z(1:n)-z(n+1:end)));
    out.nMat  = out.nMat + 2;

    mu        = z'*s/n2;
    rel_dgap  = mu/(1.0 + abs(c'*z + 0.5*z'*[ATAblkmat;-ATAblkmat] + 0.5*norm(b)^2));

    dual_e    = s - c - [ATAblkmat;-ATAblkmat];
    dual_inf  = norm(dual_e)/normc;
    
    % Print info.
    if opts.verbose == 2
        fprintf('  %5s | %3d   | %1.2e | %1.2e | \n',...
        'final', out.iters, rel_dgap, dual_inf);
    end   
    
    if (rel_dgap < opts.tol) 
        x = z(1:n) - z(n+1:end);
        out.status  = 1;
        out.message = 'Converged'; 
        break;
    elseif (out.iters > opts.maxiters)
        x = z(1:n) - z(n+1:end);
        out.status  = 0;
        out.message = 'Maximum iterations reached'; 
        break;
    end

    out.iters = out.iters + 1;
end

out.linsystems = out.linsystems + out.correctors + out.predictors;
tend = cputime;
out.totcputime = tend - tstart;

if opts.verbose >= 1
   fprintf('Total CPU time: %f \n', out.totcputime);
   fprintf('CPU time for stepsizes: %f \n', out.cpusteps);
   fprintf('Iterations: %3d \n', out.iters);
   fprintf('Predictors: %3d \n', out.predictors);
   fprintf('Correctors: %3d \n', out.correctors);
   fprintf('Total linear systems solved: %3d \n', out.linsystems);
   fprintf('Total PCG iterations: %4d \n', out.totpcgiters);
   fprintf('Total nMat: %5d \n', out.nMat);
end

end
