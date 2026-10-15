%% Dependence of the PolyOpt results on M
% This script computes the optimized step sizes for the two-scale cluster
% spectrum with fast radius r_fast = 5, for different values of M.
% It reproduces Tables 4.13-4.16 in Section 4.3.2.
thisFile = matlab.desktop.editor.getActiveFilename;
baseDir = fileparts(thisFile);
cd(baseDir);

repoDir = fullfile(baseDir,'..','..','..','Extern-Repos','RK-Opt-add_2N_fixed');
addpath(genpath(repoDir),'-begin');
clear functions

outdir = fullfile(baseDir,'results_M_dependence');
if ~exist(outdir,'dir')
    mkdir(outdir);
end


%% Parameters

basis_list = { ...
    'chebyshev', ...
    'chebinterp'};
M_list = [300 500 1000 2000 3000 4000];
kappa = 0.05;


%% 

p = 1;
s_list = 7:9;

[T,Tsum] = run_M_dependence( ...
    p,s_list,basis_list,M_list,kappa,outdir);

%% 
p = 2;
s_list = 6:8;

[T2,Tsum2] = run_M_dependence( ...
    p,s_list,basis_list,M_list,kappa,outdir);

%% 
p = 3;
s_list = 7:9;

[T3,Tsum3] = run_M_dependence( ...
    p,s_list,basis_list,M_list,kappa,outdir);

%% 
p = 4;
s_list = 8:11;

[T4,Tsum4] = run_M_dependence( ...
    p,s_list,basis_list,M_list,kappa,outdir);



%%
% Compute PolyOpt results for different values of M.
function [T,Tsum] = run_M_dependence( ...
    p,s_list,basis_list,M_list,kappa,outdir)

lambda_s = -1;
lambda_f = -1000;

r_s = kappa*abs(lambda_s);
r_f = kappa*abs(lambda_f);

rows = struct([]);

for b = 1:numel(basis_list)

    basis = basis_list{b};

    for s = s_list
        for M = M_list

            fprintf('basis=%s, p=%d, s=%d, M=%d ... ', ...
                basis,p,s,M);

            try

                theta = 2*pi*(0:M-1)/M;

                D_s = lambda_s + r_s*cos(theta) ...
                    + 1i*r_s*sin(theta);

                D_f = lambda_f + r_f*cos(theta) ...
                    + 1i*r_f*sin(theta);

                S = [D_s,D_f];

                [hstar,Rdata] = opt_poly_bisect(S,s,p,basis);

                % Evaluate the optimized polynomial on the slow and fast spectra.
                if strcmp(basis,'chebinterp')

                    % For chebinterp, Rdata contains the polynomial values
                    % at the Chebyshev interpolation nodes.
                    z = hstar*S;

                    zmin = min(real(z));
                    zmax = 0;

                    j = 0:s;
                    x = cos(pi*j/s);

                    % Map the Chebyshev nodes from [-1,1] to [zmin,0].
                    z_nodes = 0.5*((zmax+zmin) + (zmax-zmin)*x);
                    z_nodes = z_nodes(:);

                    f_nodes = Rdata(:);

                    R_slow = abs(bary_eval_cheb( ...
                        z_nodes,f_nodes,hstar*D_s));

                    R_fast = abs(bary_eval_cheb( ...
                        z_nodes,f_nodes,hstar*D_f));

                else

                    % polyval expects the coefficients in decreasing order.
                    coeffs = Rdata(:).';

                    R_slow = abs(polyval( ...
                        fliplr(coeffs),hstar*D_s));

                    R_fast = abs(polyval( ...
                        fliplr(coeffs),hstar*D_f));

                end

                max_R_slow = max(R_slow);
                max_R_fast = max(R_fast);
                max_R = max(max_R_slow,max_R_fast);

                row = struct( ...
                    'basis',string(basis), ...
                    'p',p, ...
                    's',s, ...
                    'M',M, ...
                    'hstar',hstar, ...
                    'max_R_slow',max_R_slow, ...
                    'max_R_fast',max_R_fast, ...
                    'max_R',max_R);

                rows = [rows; row]; %#ok<AGROW>

                fprintf( ...
                    'h*=%.6g, slow=%.3e, fast=%.3e, max=%.3e\n', ...
                    hstar,max_R_slow,max_R_fast,max_R);

            catch ME

                warning( ...
                    'Run failed for basis=%s, p=%d, s=%d, M=%d:\n%s', ...
                    basis,p,s,M,ME.message);

            end
            end
        end
end

T = struct2table(rows);
disp(T);


%% Summary over all values of M

summary = struct([]);

for b = 1:numel(basis_list)

    basis = string(basis_list{b});

    for s = s_list

        Tsub = T(T.basis == basis & T.s == s,:);

        if isempty(Tsub)
            continue
        end

        [~,i_min] = min(Tsub.hstar);
        [~,i_max] = max(Tsub.hstar);

        max_R_over_M = max(Tsub.max_R);

        row = struct( ...
            'basis',basis, ...
            'p',p, ...
            's',s, ...
            'h_min',Tsub.hstar(i_min), ...
            'h_max',Tsub.hstar(i_max), ...
            'max_R_slow_at_h_max',Tsub.max_R_slow(i_max), ...
            'max_R_fast_at_h_max',Tsub.max_R_fast(i_max), ...
            'max_R_over_M',max_R_over_M );

        summary = [summary; row]; %#ok<AGROW>

    end
end

Tsum = struct2table(summary);
disp(Tsum);

% Save results

outputFile = fullfile(outdir, ...
    sprintf('M_dependence_kappa%g_p%d.xlsx',kappa,p));

writetable(T,outputFile, ...
    'Sheet','Full results', ...
    'WriteMode','overwritesheet');

writetable(Tsum,outputFile, ...
    'Sheet','Summary', ...
    'WriteMode','overwritesheet');

fprintf('\nResults saved in:\n%s\n',outputFile);

end


%% Helper functions
% Evaluate the interpolation polynomial in barycentric form.
function p = bary_eval_cheb(z_nodes, f_nodes, z_eval)

    z_nodes = z_nodes(:);
    f_nodes = f_nodes(:);
    z_eval  = z_eval(:);

    N = length(z_nodes) - 1;

    w = ones(N+1,1);
    w(1)   = 0.5;
    w(end) = 0.5;
    w = w .* (-1).^(0:N)';

    p = zeros(size(z_eval));

    for k = 1:length(z_eval)
        zk = z_eval(k);

        idx = find(abs(zk - z_nodes) < 1e-14, 1);

        if ~isempty(idx)
            p(k) = f_nodes(idx);
        else
            temp = w ./ (zk - z_nodes);
            p(k) = sum(temp .* f_nodes) / sum(temp);
        end
    end

end



