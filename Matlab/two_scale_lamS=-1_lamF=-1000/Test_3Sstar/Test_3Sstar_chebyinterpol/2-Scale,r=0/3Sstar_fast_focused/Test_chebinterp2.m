% Tested configurations:
% p = 1: s = 3:5
% p = 2: s = 3:5
% p = 3: s = 4:5
% p = 4: s = 5:6
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

repoDir = fullfile(baseDir,'..','..','..','..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
addpath(genpath(repoDir),'-begin');
clear functions

%% Parameters

alpha = 1;
epsi  = 1e-3;

lam1 = -alpha;
lam2 = -1/epsi;
S = [lam1,lam2];

kappa = 0.0;
M     = 1200;
basis       = 'chebinterp';
method_type = '3Sstar';

%% 

p = 1;
s = 3;

n_grid = 1;
h_min  = 18072;
h_max  = 18072;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%% 
p = 1;
s = 4;

n_grid = 1;
h_min  = 32128;
h_max  = 32128;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 1;
s = 5;

n_grid = 1;
h_min  = 50200;
h_max  = 50200;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 2;
s = 3;

n_grid = 1;
h_min  = 2;
h_max  = 2;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 2;
s = 4;

n_grid = 1;
h_min  = 2625;
h_max  = 2625;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 2;
s = 5;

n_grid = 1;
h_min  = 12725;
h_max  = 12725;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 3;
s = 4;

n_grid = 1;
h_min  = 2.5130;
h_max  = 2.5130;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 3;
s = 5;

n_grid = 1;
h_min  = 2.5159;
h_max  = 2.5159;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 4;
s = 5;

n_grid = 20;
h_min  = 2.70 ;
h_max  = 2.7840 ;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
p = 4;
s = 6;

n_grid = 1;
h_min  = 2.7840 ;
h_max  = 2.7840 ;

num_starting_points = 300;

[rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points);

check_rk_stability(rk,h_star,S);
%%
% Optimize the stability polynomial with PolyOpt and construct a corresponding
% low-storage Runge-Kutta method.
function [rk,h_star] = run_linesearch_RK_optimization( ...
    baseDir,lam1,lam2,kappa,M,p,s,basis,method_type, ...
    n_grid,h_min,h_max,num_starting_points)

    % Construct the spectrum.
    r1 = kappa*abs(lam1);
    r2 = kappa*abs(lam2);

    theta = linspace(0, 2*pi, M);

    D1 = lam1 + r1*cos(theta) + 1i*r1*sin(theta);
    D2 = lam2 + r2*cos(theta) + 1i*r2*sin(theta);
    S  = [D1,D2];


    % Optimize the stability polynomial by line search
    [h_star,R_vals] = opt_poly_linesearch( ...
        S,s,p,basis,n_grid, ...
        'h_min',h_min, ...
        'h_max',h_max);


    % Construct the Chebyshev interpolation nodes corresponding to the
    % polynomial values returned by PolyOpt.
    z      = h_star*S;
    zmin   = min(real(z));
    zmax   = 0;
    N      = s;
    j      = 0:N;
    x_cheb = cos(pi*j/N);

    z_nodes_full = 0.5*((zmax+zmin) + ...
                        (zmax-zmin)*x_cheb);

    z_nodes_full = z_nodes_full(:);
    R_vals_full  = R_vals(:);

    % Select K = s-p interpolation nodes closest to the fast eigenvalue
    K = s-p;
    z_fast_center = h_star*lam2;

    [~,sort_idx] = sort(abs(z_nodes_full-z_fast_center));

    idx      = sort_idx(1:K);
    idx      = sort(idx);
    z_nodes  = z_nodes_full(idx);
    R_vals_K = R_vals_full(idx);

    % Evaluate the PolyOpt polynomial on the scaled spectrum
    R_eval = bary_eval_cheb( ...
        z_nodes_full,R_vals_full,h_star*S);

    mx = max(abs(R_eval));

    fprintf('h = %.6g,   max |R(h*lambda)| = %.3e\n', ...
        h_star,mx);


    % RK-Opt to construct the corresponding RK-Method
    rk = rk_opt(s,p,method_type,'acc', ...
        'interp_nodes',z_nodes, ...
        'interp_values',R_vals_K, ...
        'algorithm','interior-point', ...
        'num_starting_points',num_starting_points, ...
        'write_to_file',1, ...
        'append_time',1, ...
        'suppress_warnings',true);


    % Save the data

    dataDir = fullfile(baseDir,'polyopt_data');

    if ~exist(dataDir,'dir')
        mkdir(dataDir);
    end


    files = dir(fullfile( ...
        baseDir,sprintf('%s-%d-%d_*.txt',method_type,p,s)));

    if isempty(files)
        error(['No method file of the form %s-%d-%d_*.txt was found in:\n%s'], ...
            method_type,p,s,baseDir);
    end


    [~,idxNewest] = max([files.datenum]);
    methodName = files(idxNewest).name;

    [~,baseName,~] = fileparts(methodName);
    dataFile = fullfile(dataDir,[baseName '_plot_data.mat']);


    experiment = struct();

    experiment.method_file = methodName;
    experiment.method_type = method_type;
    experiment.basis = basis;

    experiment.lam1 = lam1;
    experiment.lam2 = lam2;
    experiment.kappa = kappa;
    experiment.M = M;
    experiment.S = S(:);

    experiment.p = p;
    experiment.s = s;

    experiment.h_min = h_min;
    experiment.h_max = h_max;
    experiment.n_grid = n_grid;
    experiment.h_star = h_star;

    experiment.z_nodes_full = z_nodes_full;
    experiment.R_vals_full = R_vals_full;


    save(dataFile,'experiment');

    fprintf('Saved plot data file:\n%s\n',dataFile);

end

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



%% Stability check
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


