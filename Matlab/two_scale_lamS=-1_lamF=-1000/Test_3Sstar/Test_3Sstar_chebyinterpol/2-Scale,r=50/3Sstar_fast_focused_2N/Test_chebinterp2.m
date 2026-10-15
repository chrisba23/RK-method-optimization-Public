% Tested configurations:
% p = 4: s = 5
% 5 runs were performed for s = 5 but none of them returned a method
thisFile = matlab.desktop.editor.getActiveFilename;
baseDir = fileparts(thisFile);
cd(baseDir);



repoDir = fullfile(baseDir,'..','..','..','..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
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
p = 4;
s = 5;

% Optimize the stability with PolyOpt
[h_star, R_vals] = opt_poly_bisect(S, s, p, 'chebinterp');
fprintf('p=%d, s=%d:  h* = %.6g', p, s, h_star);

% Construct the Chebyshev interpolation nodes corresponding to the
% polynomial values returned by PolyOpt.
z      = h_star * S;
zmin   = min(real(z));
zmax   = 0;
N      = s;
j      = 0:N;
x_cheb = cos(pi*j/N);
z_nodes_full = 0.5*((zmax+zmin) + (zmax - zmin)*x_cheb);
z_nodes_full = z_nodes_full(:);   
R_vals_full = R_vals(:);         

% Select K = s-p interpolation nodes closest to the fast eigenvalue.             
K = s-p; 
z_fast_center = h_star * lam2;

[~, sort_idx] = sort(abs(z_nodes_full - z_fast_center));
idx      = sort_idx(1:K);
idx      = sort(idx);
z_nodes  = z_nodes_full(idx);
R_vals_K = R_vals_full(idx);


% Evaluate the PolyOpt polynomial on the scaled spectrum.
R_eval = bary_eval_cheb(z_nodes_full, R_vals_full, h_star*S);
mx = max(abs(R_eval));
fprintf(',   max |R(h*lambda)| = %.3e\n', mx);


% Construct the corresponding low-storage Runge-Kutta method with RK-Opt.
rk = rk_opt(s, p, '2N', 'acc', ...
    'interp_nodes',  z_nodes, ...
    'interp_values', R_vals_K, ...
    'algorithm','interior-point', ...
    'num_starting_points', 300 , ...
    'write_to_file', 1, ...
    'append_time', 1, ...
    'suppress_warnings', true);

check_rk_stability(rk,h_star,S);

% Save the data.
dataDir = fullfile(baseDir, 'polyopt_data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

method_type = '3Sstar'; 
basis = 'chebinterp';

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

experiment.z_nodes_full = z_nodes_full;
experiment.R_vals_full = R_vals_full;

save(dataFile, 'experiment');

fprintf('Saved plot data file:\n%s\n', dataFile);

%% Helper functions

% Evaluate the interpolating polynomial in barycentric form.
function p_eval = bary_eval_cheb(z_nodes,f_nodes,z_eval)

    z_nodes = z_nodes(:);
    f_nodes = f_nodes(:);
    z_eval  = z_eval(:);

    N = length(z_nodes)-1;

    w = ones(N+1,1);
    w(1)   = 0.5;
    w(end) = 0.5;
    w = w.*(-1).^(0:N)';

    p_eval = zeros(size(z_eval));

    for k = 1:length(z_eval)

        zk = z_eval(k);

        idx = find(abs(zk-z_nodes) < 1e-14,1);

        if ~isempty(idx)
            p_eval(k) = f_nodes(idx);
        else
            temp = w./(zk-z_nodes);
            p_eval(k) = sum(temp.*f_nodes)/sum(temp);
        end

    end
end

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