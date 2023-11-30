%% Analysis of the results
% The script uses matlab2tikz to produce figures:
% https://github.com/matlab2tikz/matlab2tikz
% and some bash-command to parse the script, thus it can be run only on a
% Linux machine.

clear; clc; close all;
addpath("../../matlab2tikz/src/");

%% V-SVBM-4CHEBY4-30l1JAC
[~,logfiles] = system('ls -v1 soc1/log_vsvbm_cheby4_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_svbm_cheby4 = maketable(logfiles);

%% V-SVBM-4OPTCHEBY4-30l1JAC
[~,logfiles] = system('ls -v1 soc1/log_vsvbm_optcheby4_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_svbm_optcheby4 = maketable(logfiles);

%% V-SVBM-4OPTCHEBY1-30l1JAC
[~,logfiles] = system('ls -v1 soc1/log_vsvbm_optcheby1_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_svbm_optcheby1 = maketable(logfiles);

%% Plot figures
figure("Position",[1032 798 871 418])
hold on
loglog(tab_svbm_cheby4,"np","it","LineWidth",2,"Marker","+","DisplayName","Chebyshev 4\textsuperscript{th} kind")
loglog(tab_svbm_optcheby4,"np","it","LineWidth",2,"Marker","v","DisplayName","Optimized Chebyshev 4\textsuperscript{th} kind")
loglog(tab_svbm_optcheby1,"np","it","LineWidth",2,"Marker","o","DisplayName","Optimized Chebyshev 1\textsuperscript{st} kind")
hold off
xticks(tab_svbm_cheby4.np(~isnan(tab_svbm_cheby4.np)))
xlabel("MPI tasks")
ylabel("Iteration")
title("VBM Aggregation - $\ell_1$$-Jacobi coarse solver")
set(gca,"XScale","log")
set(gca,"YScale","log")
legend('Location','westoutside')
matlab2tikz('filename','soc1_iteration_toeplitz.tex', ...
    'width','0.45\columnwidth', ...
    'height','1.5in', ...
    'parseStrings',false)

figure("Position",[1032 798 871 418])
hold on
loglog(tab_svbm_cheby4,"np","t_solve","LineWidth",2,"Marker","+","DisplayName","Chebyshev 4\textsuperscript{th} kind")
loglog(tab_svbm_optcheby4,"np","t_solve","LineWidth",2,"Marker","v","DisplayName","Optimized Chebyshev 4\textsuperscript{th} kind")
loglog(tab_svbm_optcheby1,"np","t_solve","LineWidth",2,"Marker","o","DisplayName","Optimized Chebyshev 1\textsuperscript{st} kind")
hold off
xticks(tab_svbm_cheby4.np(~isnan(tab_svbm_cheby4.np)))
xlabel("MPI tasks")
ylabel("Solve Time (s)")
title("VBM Aggregation - $\ell_1$--Jacobi coarse solver")
set(gca,"XScale","log")
set(gca,"YScale","log")
legend('Location','westoutside')
matlab2tikz('filename','soc1_solvetime_toeplitz.tex', ...
    'width','0.45\columnwidth', ...
    'height','1.5in', ...
    'parseStrings',false)

figure("Position",[1032 798 871 418])
hold on
loglog(tab_svbm_cheby4,"np","t_it","LineWidth",2,"Marker","+","DisplayName","Chebyshev 4\textsuperscript{th} kind")
loglog(tab_svbm_optcheby4,"np","t_it","LineWidth",2,"Marker","v","DisplayName","Optimized Chebyshev 4\textsuperscript{th} kind")
loglog(tab_svbm_optcheby1,"np","t_it","LineWidth",2,"Marker","o","DisplayName","Optimized Chebyshev 1\textsuperscript{st} kind")
hold off
xticks(tab_svbm_cheby4.np(~isnan(tab_svbm_cheby4.np)))
xlabel("MPI tasks")
ylabel("Time per Iteration (s)")
title("VBM Aggregation - $\ell_1$-Jacobi coarse solver")
set(gca,"XScale","log")
set(gca,"YScale","log")
legend('Location','westoutside')
matlab2tikz('filename','soc1_titeration_toeplitz.tex', ...
    'width','0.45\columnwidth', ...
    'height','1.5in', ...
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