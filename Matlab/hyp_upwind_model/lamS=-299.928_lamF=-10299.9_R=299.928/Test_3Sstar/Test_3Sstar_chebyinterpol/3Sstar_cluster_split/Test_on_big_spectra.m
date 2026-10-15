% Tested configurations:
% p = 3: s = 9:14
% Each configuration was run once in RK-Opt.

thisFile = matlab.desktop.editor.getActiveFilename;
baseDir = fileparts(thisFile);
cd(baseDir);

repoDir = fullfile(baseDir,'..','..','..','..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
addpath(genpath(repoDir),'-begin');
clear functions
%% Parameters

% Define the centers and radius of the spectrum.
lamS = -299.928;
lamF = -10299.9;
R    = 299.928;

% Construct the spectrum.
M = 800;                        
theta = linspace(0, 2*pi, M);

Dslow = lamS + R*cos(theta) + 1i*R*sin(theta);
Dfast = lamF + R*cos(theta) + 1i*R*sin(theta);

S = [Dslow, Dfast];

%% Optimization
p = 3;
s = 13;

% Optimize the stability with PolyOpt
[h_star, R_vals] = opt_poly_bisect(S, s, p, 'chebinterp');

fprintf('Case 1: p=%d, s=%d, h* = %.8g\n', p, s, h_star);

% Construct the Chebyshev interpolation nodes corresponding to the
% polynomial values returned by PolyOpt.
z      = h_star * S;
zmin   = min(real(z));
zmax   = 0;
N      = s;
j      = 0:N;
x_cheb = cos(pi*j/N);
z_nodes_full = 0.5*((zmax+zmin) + (zmax - zmin)*x_cheb);
z_nodes_full = z_nodes_full(:);   % (s+1) x 1
R_vals_full = R_vals(:);         


% Split K = s-p interpolation nodes between the slow and fast cluster.  
K = s - p;
K1 = floor(K/2);   
K2 = K - K1;       

j1 = (1:K1)';
x1 = cos(pi*j1/(K1+1));

zmin1 = h_star * (real(lamS) - R);
zmax1 = h_star * (real(lamS) + R);

z_nodes1 = 0.5*((zmax1+zmin1) + (zmax1-zmin1)*x1);


j2 = (1:K2)';
x2 = cos(pi*j2/(K2+1));

zmin2 = h_star * (real(lamF) - R);
zmax2 = h_star * (real(lamF) + R);

z_nodes2 = 0.5*((zmax2+zmin2) + (zmax2-zmin2)*x2);


z_nodes = sort([z_nodes1; z_nodes2]);
R_vals_K = bary_eval_cheb(z_nodes_full, R_vals_full, z_nodes);

z_nodes_full = z_nodes_full(:);
R_vals_full = R_vals(:);

% Evaluate the PolyOpt polynomial on the scaled spectrum.
R_eval = bary_eval_cheb(z_nodes_full, R_vals_full, h_star*S);
mx = max(abs(R_eval));
fprintf(',   max |R(h*lambda)| = %.3e\n', mx);

% Construct the corresponding low-storage Runge-Kutta method with RK-Opt.
rk = rk_opt(s, p, '3Sstar', 'acc', ...
    'interp_nodes',  z_nodes, ...
    'interp_values', R_vals_K, ...
    'algorithm','interior-point', ... 
    'num_starting_points', 200, ...
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