%% Analysis of the results
% The script uses matlab2tikz to produce figures:
% https://github.com/matlab2tikz/matlab2tikz
% and some bash-command to parse the script, thus it can be run only on a
% Linux machine.

clear; clc; close all;
addpath("../../matlab2tikz/src/");

TSOLVE = zeros(3,2);
TSOLVE2 = zeros(3,2);
TSOLVE3 = zeros(3,2);

%% V-SVBM-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 soc1/log_POLY_LOTTES_svbm*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_svbm_cheby4 = maketable(logfiles);

TSOLVE(1,:) = tab_svbm_cheby4.t_solve(:);

%% V-SVBM-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 soc1/log_POLY_LOTTES_BETA_svbm*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_svbm_optcheby4 = maketable(logfiles);

TSOLVE(2,:) = tab_svbm_optcheby4.t_solve(:);

%% V-SVBM-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 soc1/log_POLY_NEW_svbm*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_svbm_otcheby1 = maketable(logfiles);

TSOLVE(3,:) = tab_svbm_otcheby1.t_solve(:);

%% V-SMATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 match/log_POLY_LOTTES_smatch*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_cheby4 = maketable(logfiles);

TSOLVE2(1,:) = tab_match_cheby4.t_solve(:);

%% V-SMATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 match/log_POLY_LOTTES_BETA_smatch*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby4 = maketable(logfiles);

TSOLVE2(2,:) = tab_match_optcheby4.t_solve(:);

%% V-SMATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 match/log_POLY_NEW_smatch*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_otcheby1 = maketable(logfiles);

TSOLVE2(3,:) = tab_match_otcheby1.t_solve(:);

%% V-KMATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 match/log_POLY_LOTTES_kmatch*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_kmatch_cheby4 = maketable(logfiles);

TSOLVE3(1,:) = tab_kmatch_cheby4.t_solve(:);

%% V-KMATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 match/log_POLY_LOTTES_BETA_kmatch*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_kmatch_optcheby4 = maketable(logfiles);

TSOLVE3(2,:) = tab_kmatch_optcheby4.t_solve(:);

%% V-KMATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 match/log_POLY_NEW_kmatch*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_kmatch_otcheby1 = maketable(logfiles);

TSOLVE3(3,:) = tab_kmatch_otcheby1.t_solve(:);

%% Plot
X = [1,2,3];
figure("Position",[197 532 1397 387])
subplot(1,3,1)
b = bar(X,TSOLVE);
b(1).CData = repmat([0 104 181]/255,3,1); 
b(2).CData = repmat([118 185 0]/255,3,1); 
legend({'i7-8750H','GeForce GTX 1050'},"Location","northeast")
title("VBM Aggregation")
xticklabels({'Chebyshev 4','Opt. Chebyshev 4','Opt. Chebyshev 1'});

subplot(1,3,2)
b = bar(X,TSOLVE2);
b(1).CData = repmat([0 104 181]/255,3,1); 
b(2).CData = repmat([118 185 0]/255,3,1); 
title("(Smoothed) Matching Aggregation")
xticklabels({'Chebyshev 4','Opt. Chebyshev 4','Opt. Chebyshev 1'});

subplot(1,3,3)
b = bar(X,TSOLVE3);
b(1).CData = repmat([0 104 181]/255,3,1); 
b(2).CData = repmat([118 185 0]/255,3,1); 
title("(Unsmoothed) Matching Aggregation")
xticklabels({'Chebyshev 4','Opt. Chebyshev 4','Opt. Chebyshev 1'});

matlab2tikz('filename','x580gd_gpu.tex','width','\columnwidth', ...
    'parseStrings',false)

%% Function to analyze logs
function tab = maketable(logfiles)
%MAKETABLE produces tables with the result of a given run
n = length(logfiles);

for i = 1:n
    file = logfiles{i};

    [~,out] = system(sprintf('grep "Iterations to convergence          :" %s',file));
    try
        iter(i,1) = sscanf(out,'Iterations to convergence          : %d');
    catch
        iter(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Time to build hierarchy            :" %s',file));
    try
        thierarchy(i,1) = sscanf(out,'Time to build hierarchy            : %e');
    catch
        thierarchy(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Time to build smoothers            :" %s',file));
    try
        tsmoothers(i,1) = sscanf(out,'Time to build smoothers            : %e');
    catch
        tsmoothers(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Total time for preconditioner      :" %s',file));
    try
        tbuild(i,1) = sscanf(out,'Total time for preconditioner      : %e');
    catch
        tbuild(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Time to solve system               :" %s',file));
    try
        tsolve(i,1) = sscanf(out,'Time to solve system               : %e');
    catch
        tsolve(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Time per iteration                 :" %s',file));
    try
        titer(i,1) = sscanf(out,'Time per iteration                 : %e');
    catch
        titer(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Total time                         :" %s',file));
    try
        tottime(i,1) = sscanf(out,'Total time                         : %e');
    catch
        tottime(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Linear system size                 : " %s',file));
    try
        sysize(i,1) = sscanf(out,'Linear system size                 : %e');
    catch
        sysize(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "   Number of levels   :            " %s',file));
    try
        levels(i,1) = sscanf(out,'   Number of levels   : %d');
    catch
        levels(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "   Operator complexity:    " %s',file));
    try
        opc(i,1) = sscanf(out,'   Operator complexity: %f');
    catch
        opc(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "   Average coarsening :    " %s',file));
    try
        avgcoars(i,1) = sscanf(out,'   Average coarsening : %f');
    catch
        avgcoars(i,1) = NaN;
    end

    [~,out] = system(sprintf('grep "Computed solution on " %s',file));
    try
        np(i,1) = sscanf(out,'Computed solution on %d');
    catch
        np(i,1) = NaN;
    end
end

tab = table(np,sysize,thierarchy,tsmoothers,tbuild,tsolve,iter,titer, ...
    opc,avgcoars,levels,'VariableNames',{'np','n','t_hier','t_smooth', ...
    't_build','t_solve','it','t_it','opc','avg_coars','lev'});

end