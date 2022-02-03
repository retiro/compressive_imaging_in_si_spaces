%% Compressive Imaging in Shift-Invariant Spaces
%
% Implementation of
% T. Vlasic and D. Sersic, "Single-pixel compressive imaging in shift-
% invariant spaces via exact wavelet frames," 2021, arXiv:2106.00404.
%
% MATLAB script
% 
% Copyright (c) 2021, Tin Vlasic

%%
close all
clearvars
clc
rng('default')

%% Definition and initialization of variables
p_in = 3; % B-spline degree of the input signal
N = 512; % number of B0-spline sampling functions in acquisition process 
         % in 1D (exponent with base 2)

% B-spline direct and indirect transform variables
boundExtension = 'symmetric'; % boundary extension type
precision = 1e-4; % precision of B-spline representation
upSamp = 1; % upsampling factor when interpolating a reconstructed image
            % reconstructructed image size (Nxm)x(Nxm)

% Configure measurements
meas_ratio = 0.25; % measurement ratio
measMtxType = 'wht'; % 'dct', 'fft' or 'wht'
subsampling = 'random'; % 'random' or 'multilevel'

% Sparsity
sparsityDomain = 'wavelets'; % 'dct' or 'wavelets'

% Wavelet type and variables
scaleLevels = 3; % define the number of scale levels
wave_type = 'biorthogonal'; % 'biorthogonal' or 'orthogonal'
wav_bior = 'bior2.2'; % wavelet function, I.J - I should be p_in+1
                      % waveinfo('bior')
                      % JPEG - bior2.2 (CDF5/3) and bior4.4 (CDF9/7)
                      % I = 1, J = 1, 3, 5
                      % I = 2, J = 2, 4, 6, 8
                      % I = 3, J = 1, 3, 5, 7, 9

% Load test image
imageFolderName = ['.', filesep, 'standard_test_images', filesep];
imName = 'cameraman.tif';
% imName = 'lena_gray_512.tif';
% imName = 'peppers_gray.tif';
% imName = 'pirate.tif';
% imName = 'barbara.png';
% imName = 'boat.tif';

% Fix the seed
n_rng = 2; % in paper we use 2
rng(n_rng)
                      
%% global variables (do not edit)
n_tilde = floor(p_in/2);
N_overlap_sig = floor((p_in + 1)/2); % number of overlaping basis functions
M = floor(meas_ratio*N^2); % number of measurements, ratio = M/N^2

fprintf('Measurement resolution: %d x %d\n', N, N);
fprintf('Measurement ratio: m=%d/%d=%1.2f\n', M, N^2, meas_ratio);
if strcmp(sparsityDomain, 'wavelets')
    fprintf('Sparsity matrix: %s %s\n', wave_type, sparsityDomain);
else
    fprintf('Sparsity matrix: %s\n', sparsityDomain);
end

% maximal scale level for biorthogonal wavelets
if contains(wave_type, 'orthogonal') && strcmp(sparsityDomain, 'wavelets')
    L = wmaxlev(N, wav_bior);
    if L < scaleLevels
        scaleLevels = L;
    end
    fprintf('Wavelet function: %s\n', wav_bior);
    fprintf('Scale/maxScale: %d/%d\n', scaleLevels, L);
end
fprintf('B-spline generator phi(t) of order: %d\n', p_in);
fprintf('B-spline sampling kernel zeta(t) of order: 0\n');

%% Calculate cross-correlation sequence
r_sa = crossCorr(p_in, 0);

%% Preprocess image
fprintf('Input image: %s\n', imName);
Im = im2double(imread([imageFolderName, imName]));
try Im = rgb2gray(Im); catch; end
if size(Im,3)>1; Im = Im(:,:,1); end
fprintf('Input image size: %d x %d pixels\n', size(Im,1), size(Im,2));

fprintf('Upsampling factor: %d\n', upSamp);
Im = imresize(Im, [N*upSamp, N*upSamp], 'cubic');
fprintf('Rescaled input image size: (%dx%d)x(%dx%d)= %d x %d pixels\n', N, upSamp, N, upSamp, size(Im,1), size(Im,2));

% plot image
figure('Renderer', 'painters', 'Position', [100 100 1000 400]);
subplot(121); imagesc(Im); colormap gray; title('Ground Truth');
cl = caxis; % save c-axis scale

%% Simulate CS measurement process
if strcmp(subsampling, 'random')
    % uniform random subsampling
    picks = sort(randsample(1:N*N, M, false))'; % subsampled rows
    picks(1) = 1; % measure DC
    perm = randperm(N*N)'; % global randomizer
    [~, permReturn] = sort(perm); % return sequence of permutation
elseif strcmp(subsampling, 'multilevel')
    % multilevel subsampling
    mLevPicks = multiLevSubsample( N*N, scaleLevels, M); % subsampled rows
    picks = [mLevPicks{scaleLevels+2:-1:1}]';
    picks(1) = 1; % measure DC
else
    error('Wrong subsampling method. (Choose between ''random'' and ''multilevel'').');
end

% DMD resolution
funSum = @(x) sum(sum(x.data))/upSamp/upSamp; 
Im_MeasRes = blockproc(Im, [upSamp, upSamp], funSum);
Im_MeasRes = Im_MeasRes(:);

% measurement procedure
if strcmp(measMtxType, 'dct')
    % Phi - random subsampling, discrete time transform, global randomizer
    Im_MeasRes_perm = Im_MeasRes(perm); % randomize sample locations
    y_all = dct(Im_MeasRes_perm); % full NxN measurements
    y = y_all(picks); % reduced set of measurements
    if strcmp(subsampling, 'multilevel')
        warning('DCT measurement matrix implementation uses only random subsampling method.');
    end
elseif strcmp(measMtxType, 'fft')
    % Phi - random subsampling, fast Fourier transform, global randomizer
    Im_MeasRes_perm = Im_MeasRes(perm); % randomize sample locations
    y_all = fft(Im_MeasRes_perm)./N; % full NxN measurements
    y = y_all(picks); % reduced set of measurements
    if strcmp(subsampling, 'multilevel')
        warning('FFT measurement matrix implementation uses only random subsampling method.');
    end
elseif strcmp(measMtxType, 'wht')
    if strcmp(subsampling, 'random')
        % Phi - random subsampling, Walsh-Hadamard transform, global randomizer
        Im_MeasRes_perm = Im_MeasRes(perm); % randomize sample locations
        y_all = fastWht(Im_MeasRes_perm')'; % full NxN measurements
        y = y_all(picks); % reduced set of measurements
    elseif strcmp(subsampling, 'multilevel')
        y_all = fastWht(Im_MeasRes')'; % full NxN measurements
        y = y_all(picks); % reduced set of measurements
    else
        error('Wrong subsampling method. (Choose between ''random'' and ''multilevel'').');
    end
else
    error('Wrong measurement matrix type or not selected. (Choose between ''dct'', ''fft'' or ''wht'').');
end

%% solve l1-minimization problem
% min tau*||x||_1 + 0.5*||Ax - y||^2,
fprintf('\nSolver details\n');

if strcmp(sparsityDomain, 'wavelets')
    % generate bookkeeping matrix S for wavelet decomposition and reconstruction
    [x_tmp, S] = wavedec2(ones(N+2*N_overlap_sig, N+2*N_overlap_sig), scaleLevels, wav_bior);
    x_len = numel(x_tmp);
    clear x_tmp;
    % use matrix-free interior-point method (mfIPM)
    if strcmp(measMtxType, 'dct')
        A = A_operator( @(x) mfSensOperator_randSubsampDctGlobalRand_dwt(x,S,1,picks,perm,permReturn,r_sa,wav_bior) ,...
                        @(x) mfSensOperator_randSubsampDctGlobalRand_dwt(x,S,2,picks,perm,permReturn,r_sa,wav_bior) );
    elseif strcmp(measMtxType, 'fft')
        A = A_operator( @(x) mfSensOperator_randSubsampFftGlobalRand_dwt(x,S,1,picks,perm,permReturn,r_sa,wav_bior) ,...
                        @(x) mfSensOperator_randSubsampFftGlobalRand_dwt(x,S,2,picks,perm,permReturn,r_sa,wav_bior) );
    elseif strcmp(measMtxType, 'wht')
        if strcmp(subsampling, 'random')
            A = A_operator( @(x) mfSensOperator_randSubsampWhtGlobalRand_dwt(x,S,1,picks,perm,permReturn,r_sa,wav_bior) ,...
                            @(x) mfSensOperator_randSubsampWhtGlobalRand_dwt(x,S,2,picks,perm,permReturn,r_sa,wav_bior) );
        elseif strcmp(subsampling, 'multilevel')
            A = A_operator( @(x) mfSensOperator_multilevSubsampWht_dwt(x,S,1,picks,N*N,r_sa,wav_bior) ,...
                            @(x) mfSensOperator_multilevSubsampWht_dwt(x,S,2,picks,N*N,r_sa,wav_bior) );
        else
            error('Wrong subsampling method. (Choose between ''random'' and ''multilevel'').');
        end
    else
        error('Wrong measurement matrix type or not selected. (Choose between ''fft'' or ''wht'').');
    end
elseif strcmp(sparsityDomain, 'dct')
    x_len = (N+2*N_overlap_sig) .* (N+2*N_overlap_sig);
    % use matrix-free interior-point method (mfIPM)
    A = A_operator( @(x) mfSensOperator_randSubsampWhtGlobalRand_dct(x,1,picks,perm,permReturn,r_sa) ,...
                    @(x) mfSensOperator_randSubsampWhtGlobalRand_dct(x,2,picks,perm,permReturn,r_sa) );
    warning('DCT sparsity matrix uses only structurally random WHT measurement matrix.');
else
    error('Wrong sparsity domain or not selected. (Choose between ''dct'' or ''wavelets'').');
end
                      
% Options.
opts.tol       = 1.0e-8; % in the paper we use 1.0e-8
opts.maxiters  = 100;
opts.tolpcg    = 1.0e-2; % in the paper we use 1.0e-2
opts.mxiterpcg = 200;
opts.verbose   = 1; % 0 for no printing info, 1 for basic printing info
                    % and 2 for detailed printing info.
% ***************************************************************************************
% If the measurements are noiseless please use opts.tol >= 1.e-8. In case that the
% measurements are contaminated with noise then please use opts.tol >= 1.0e-6. Regarding
% the accuracy of the PCG method opts.tolpcg=1.0e-1 to 1.0e-4 should be enough.
% ***************************************************************************************

% Set LASSO regularization parameter
lambdaFolderName = ['.', filesep, 'results_lambda_', measMtxType, '_tol-8-2', filesep];
imName_core = split(imName, '.');
lambdaFileName = sprintf('%s_%s_sLev%d_N%d_upSamp%d.mat', char(imName_core(1)), wav_bior, scaleLevels, N, upSamp);
try
    load([lambdaFolderName lambdaFileName]);
    lambda = C.val(p_in+1, C.measRatios==meas_ratio);
catch
    warning('Parameter lambda not found in external results. Set lambda to 1e-3.');
    lambda = 1e-3;
end
a = norm(A'*y, Inf);
tau = lambda*a; % regularization parameter.
fprintf('lambda = %e\n\n', lambda);

% Call the Matrix-free IPM solver.
[x, out] = mfipm(x_len, A, y, tau, opts);

%%  Recover SI expansion coefficients a_0[n]
if strcmp(sparsityDomain, 'wavelets')
    d_rec = waverec2(x, S, wav_bior); % recover B-spline expansion coefficients
elseif strcmp(sparsityDomain, 'dct')
    x = reshape(x, [sqrt(x_len), sqrt(x_len)]);
    d_rec = idct2(x);
else
    error('Wrong sparsity domain or not selected. (Choose between ''dct'' or ''wavelets'').');
end

%% Reconstruct the image from the obtained coefficients
% Interpolate the image from its B-spline coefficients a_0
yIm = blockproc(d_rec, [size(d_rec,1), 1], @(x) bSplineIndirectTransform(x.data, p_in, upSamp));
Im_rec = blockproc(yIm', [size(d_rec,2), 1], @(x) bSplineIndirectTransform(x.data, p_in, upSamp));
Im_rec = abs(Im_rec');

Im_psnr = psnr(Im_rec, Im);
Im_ssim = ssim(Im_rec, Im);
fprintf('\n  PSNR  |  SSIM\n');
fprintf('----------------\n');
fprintf('%2.4f | %1.4f\n', Im_psnr, Im_ssim);

% plot reconstruction
subplot(1,2,2);
imagesc(Im_rec); colormap(gray);
title(['rec. basis: B', num2str(p_in), '-spline, upsamp. factor: ', num2str(upSamp),', psnr: ', num2str(Im_psnr,4), 'dB, ssim: ', num2str(Im_ssim,4)]);
caxis(cl);
