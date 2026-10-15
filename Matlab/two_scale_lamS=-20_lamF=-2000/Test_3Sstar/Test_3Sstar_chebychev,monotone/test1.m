% Tested configurations:
% p = 3: s = 17:19
% Each configuration was run once in RK-Opt.
thisFile = matlab.desktop.editor.getActiveFilename;
baseDir = fileparts(thisFile);
cd(baseDir);

repoDir = fullfile(baseDir,'..','..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
addpath(genpath(repoDir),'-begin');
clear functions

%% Parameters


A = readmatrix('matrixN20x2.csv');

% Define the centers and radius of the spectrum.
lamS = -20;
lamF = -2000;

% Extract the spectrum of the matrix
S = eig(A).';
M = numel(S);
           
%% Optimization
p = 3;
s = 17;

% Optimize the stability with PolyOpt
[h_star, a] = opt_poly_bisect(S, s, p, 'chebyshev');

mx = max(abs(polyval(flipud(a(:)), h_star*S)));
fprintf('p=%d, s=%d:  h* = %.6g,   max |R(h*λ)| = %.3e\n', p, s, h_star, mx);

% Use RK-Opt with the polynomial obtained from PolyOpt
 function rk = rk_from_poly(poly_coeff, p, s)
    m = s-p;                      % Number of polynomial coefficients to fix
    poly_ind = (s-m+1):s;        
    poly_val = poly_coeff(poly_ind + 1);  % Corresponding coefficient values

    
     rk = rk_opt(s, p, '3Sstar', 'acc', ...
         'poly_coeff_ind', poly_ind, ...
         'poly_coeff_val', poly_val, ...
         'c_monotone', true, ...
         'algorithm','interior-point', ...
         'num_starting_points',200, ...
         'write_to_file', 1, ...
         'append_time', 1, ...
         'suppress_warnings', true);
 end

rk = rk_from_poly(a,p,s);
check_rk_stability(rk,h_star,S); 


% Save the data.

dataDir = fullfile(baseDir, 'polyopt_data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

method_type = '3Sstar';
basis = 'chebyshev';

files = dir(fullfile(baseDir, sprintf('%s-%d-%d_*.txt', method_type, p, s)));


[~, idxNewest] = max([files.datenum]);
methodName = files(idxNewest).name;

[~, baseName, ~] = fileparts(methodName);
dataFile = fullfile(dataDir, [baseName '_plot_data.mat']);

experiment = struct();

experiment.method_file = methodName;
experiment.method_type = method_type;
experiment.basis = basis;

experiment.h_star = h_star;
experiment.M = M;
experiment.S = S(:);

experiment.a = a(:);

save(dataFile, 'experiment');

fprintf('Saved plot data file:\n%s\n', dataFile);

%% Helper functions

% Evaluate the RK stability function in Butcher form on the scaled spectrum.
function max_abs_R = check_rk_stability(rk,h,S)
    A_rk = rk.A;
    b = rk.b(:);

    s = numel(b);
    e = ones(s,1);
    I = eye(s);

    z = h*S(:);
    R = zeros(size(z));

    for j = 1:numel(z)
        Y = (I-z(j)*A_rk)\e;
        R(j) = 1+z(j)*(b.'*Y);
    end

    max_abs_R = max(abs(R));

    fprintf('h = %.6g,   max |R_RK(h*lambda)| = %.6e\n', ...
        h,max_abs_R);
end