%% Additional grid-based line search
% This script performs the grid-based line search for the two-scale cluster
% spectrum with fast radius r_fast = 5 and reproduces the results for
% Table 4.17
thisFile = matlab.desktop.editor.getActiveFilename;
baseDir = fileparts(thisFile);
cd(baseDir);

repoDir = fullfile(baseDir,'..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
addpath(genpath(repoDir),'-begin');
clear functions

outdir = fullfile(baseDir,'results_linesearch');

if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% Parameters

kappa = 0.005;
n_grid = 50;

%%
p = 2;
s_list = 8;
basis = 'chebinterp';

M = 7000;

h_edges = 1:5:31;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];

Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);
%%
p = 2;
s_list = 9;
basis = 'chebinterp';

M = 7000;

h_edges = 1:5:31;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);
%%
p = 2;
s_list = 10;
basis = 'chebinterp';

M = 7000;

h_edges = 1:5:31;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);
%%
p = 2;
s_list = 11;
basis = 'chebinterp';

M = 7000;

h_edges = 1:5:31;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);
%%
p = 2;
s_list = 12;
basis = 'chebinterp';

M = 7000;
n_grid1 = 100;
h_edges = 26:5:36;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid1,outdir);
%%
p = 4;
s_list = 8;
basis = 'chebinterp';

M = 4000;
h_edges = 0:1:3;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);
%%
p = 4;
s_list = 9;
basis = 'chebinterp';

M = 4000;

h_edges = 1:1:3;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);
%%
p = 4;
s_list = 10;
basis = 'chebinterp';

M = 4000;

h_edges = 1:1:3;
intervals = [h_edges(1:end-1)' h_edges(2:end)'];


Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir);


%%
% Run the grid-based line search on each step-size interval.
function Tline = run_linesearch_intervals( ...
    p,s_list,basis,M,kappa,intervals,n_grid,outdir)

lambda_s = -1;
lambda_f = -1000;

r_s = kappa*abs(lambda_s);
r_f = kappa*abs(lambda_f);

theta = linspace(0, 2*pi, M);

D_s = lambda_s + r_s*cos(theta) + 1i*r_s*sin(theta);

D_f = lambda_f + r_f*cos(theta) + 1i*r_f*sin(theta);

S = [D_s,D_f];

rows = struct([]);

for s = s_list

    fprintf('\nTesting basis=%s, p=%d, s=%d, M=%d, kappa=%g\n', ...
        basis,p,s,M,kappa);

    for j = 1:size(intervals,1)

        h_min = intervals(j,1);
        h_max = intervals(j,2);

        fprintf('interval [%.4g, %.4g] ... ',h_min,h_max);

        try

            [h_found,~] = opt_poly_linesearch( ...
                S,s,p,basis,n_grid, ...
                'h_min',h_min, ...
                'h_max',h_max);

        catch ME

            warning( ...
                'Line search failed for p=%d, s=%d, interval=[%g,%g]:\n%s', ...
                p,s,h_min,h_max,ME.message);

            h_found = NaN;

        end

        row = struct( ...
            'basis',string(basis), ...
            'p',p, ...
            's',s, ...
            'h_min',h_min, ...
            'h_max',h_max, ...
            'n_grid',n_grid, ...
            'h_found',h_found);

        rows = [rows; row]; %#ok<AGROW>

        if isnan(h_found)
            fprintf('failed\n');
        elseif h_found > 0
            fprintf('largest feasible h = %.12g\n',h_found);
        else
            fprintf('none\n');
        end

    end
end

Tline = struct2table(rows);

disp(Tline);

outputFile = fullfile(outdir, ...
    sprintf('linesearch_%s_kappa%g_p%d_s%d.xlsx', ...
    basis,kappa,p,s_list));

writetable(Tline,outputFile, ...
    'Sheet','Intervals', ...
    'WriteMode','overwritesheet');

fprintf('\nResults saved in:\n%s\n',outputFile);

end





