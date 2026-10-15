% Tested configurations:
% p = 1: s = 4:9
% p = 2: s = 7:9
% p = 3: s = 7:11
% p = 4: s = 7:10
% Each configuration was initially run once. If no suitable method was
% obtained, up to four additional RK-Opt runs were performed.
% The stage search was stopped when larger tested stage numbers did not
% further improve the step size per right-hand side evaluation or when
% no suitable method was obtained after five runs.
%
% Files from configurations classified as unsuccessful after five runs
% are stored in "unsuccessful_runs".
thisFile = matlab.desktop.editor.getActiveFilename;
baseDir = fileparts(thisFile);
cd(baseDir);

repoDir = fullfile(baseDir,'..','..','..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
addpath(genpath(repoDir),'-begin');
clear functions
%% Parameters

alpha = 1;
epsi  = 1e-3;

lam1 = -alpha;
lam2 = -1/epsi;

% Construct the spectrum.
kappa = 0.05;           
r1 = kappa*abs(lam1);
r2 = kappa*abs(lam2);

M = 300;                
theta = linspace(0, 2*pi, M);

D1 = lam1 + r1*cos(theta) + 1i*r1*sin(theta);
D2 = lam2 + r2*cos(theta) + 1i*r2*sin(theta);
S  = [D1, D2];

%% Optimization
p = 3;
s = 9;

% Optimize the stability with PolyOpt
[h_star, a] = opt_poly_bisect(S, s, p, 'chebyshev');

mx = max(abs(polyval(fliplr(a), h_star*S)));
fprintf('p=%d, s=%d:  h* = %.6g,   max |R(h*λ)| = %.3e\n', p, s, h_star, mx);


% Construct the corresponding low-storage Runge-Kutta method with RK-Opt.
function rk = rk_from_poly(poly_coeff, p, s)
    m = s-p;              % Number of polynomial coefficients to fix
    poly_ind = (s-m+1):s;
    poly_val = poly_coeff(poly_ind + 1); % Corresponding coefficient values

    rk = rk_opt(s, p, '2N', 'acc', ...
        'poly_coeff_ind', poly_ind, ... 
        'poly_coeff_val', poly_val, ...
        'algorithm', 'interior-point', ...
        'num_starting_points', 20, ...
        'write_to_file', 1, ...
        'suppress_warnings', true);

end

rk = rk_from_poly(a,p,s);
check_rk_stability(rk,h_star,S); 


% Save the data.
dataDir = fullfile(baseDir, 'polyopt_data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

method_type = '2N';
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