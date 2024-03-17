%% Analysis of the results
% The script uses matlab2tikz to produce figures:
% https://github.com/matlab2tikz/matlab2tikz
% and some bash-command to parse the script, thus it can be run only on a
% Linux machine.

clear; clc; close all;
addpath("../../../../matlab2tikz/src/");

TSOLVE1 = zeros(4,7);    % 4 Iteration/degree 
TSOLVE2 = zeros(4,7);   % 6 Iteration/degree
TSOLVE3 = zeros(4,7);   % 8 Iteration/degree
TSOLVE4 = zeros(4,7);   % 10 Iteration/degree
TSOLVE5 = zeros(4,7);   % 12 Iteration/degree
TITER1 = zeros(4,7);    % 4 Iteration/degree 
TITER2 = zeros(4,7);   % 6 Iteration/degree
TITER3 = zeros(4,7);   % 8 Iteration/degree
TITER4 = zeros(4,7);   % 10 Iteration/degree
TITER5 = zeros(4,7);   % 12 Iteration/degree

%% V-MATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 4/log_cheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_cheby4 = maketable(logfiles);

TSOLVE1(1,:) = tab_match_cheby4.t_solve(:);
TITER1(1,:) = tab_match_cheby4.t_it(:);
NGPUS = tab_match_cheby4.np(:);
%% V-MATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 4/log_optcheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby4 = maketable(logfiles);

TSOLVE1(2,:) = tab_match_optcheby4.t_solve(:);
TITER1(2,:) = tab_match_optcheby4.t_it(:);
%% V-MATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 4/log_optcheby1_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby1 = maketable(logfiles);

TSOLVE1(3,:) = tab_match_optcheby1.t_solve(:);
TITER1(3,:) = tab_match_optcheby1.t_it(:);
%% V-MATCH-4L1-JAC-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 4/log_l1jacobi_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_l1jac = maketable(logfiles);

TSOLVE1(4,:) = tab_match_l1jac.t_solve(:);
TITER1(4,:) = tab_match_l1jac.t_it(:);
%% Plot analysis
figure(1)
subplot(3,2,1)
loglog(NGPUS,TSOLVE1(1,:), ...
    NGPUS,TSOLVE1(2,:), ...
    NGPUS,TSOLVE1(3,:), ...
    NGPUS,TSOLVE1(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,1)
loglog(NGPUS,TITER1(1,:), ...
    NGPUS,TITER1(2,:), ...
    NGPUS,TITER1(3,:), ...
    NGPUS,TITER1(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')

%% V-MATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 6/log_cheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_cheby4 = maketable(logfiles);

TSOLVE2(1,:) = tab_match_cheby4.t_solve(:);
TITER2(1,:) = tab_match_cheby4.t_it(:);
NGPUS = tab_match_cheby4.np(:);
%% V-MATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 6/log_optcheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby4 = maketable(logfiles);

TSOLVE2(2,:) = tab_match_optcheby4.t_solve(:);
TITER2(2,:) = tab_match_optcheby4.t_it(:);
%% V-MATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 6/log_optcheby1_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby1 = maketable(logfiles);

TSOLVE2(3,:) = tab_match_optcheby1.t_solve(:);
TITER2(3,:) = tab_match_optcheby1.t_it(:);
%% V-MATCH-4L1-JAC-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 6/log_l1jacobi_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_l1jac = maketable(logfiles);

TSOLVE2(4,:) = tab_match_l1jac.t_solve(:);
TITER2(4,:) = tab_match_l1jac.t_it(:);
%% Plot analysis
figure(1)
subplot(3,2,2)
loglog(NGPUS,TSOLVE2(1,:), ...
    NGPUS,TSOLVE2(2,:), ...
    NGPUS,TSOLVE2(3,:), ...
    NGPUS,TSOLVE2(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
figure(2)
subplot(3,2,2)
loglog(NGPUS,TITER2(1,:), ...
    NGPUS,TITER2(2,:), ...
    NGPUS,TITER2(3,:), ...
    NGPUS,TITER2(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on

%% V-MATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 8/log_cheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_cheby4 = maketable(logfiles);

TSOLVE3(1,:) = tab_match_cheby4.t_solve(:);
TITER3(1,:) = tab_match_cheby4.t_it(:);
NGPUS = tab_match_cheby4.np(:);
%% V-MATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 8/log_optcheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby4 = maketable(logfiles);

TSOLVE3(2,:) = tab_match_optcheby4.t_solve(:);
TITER3(2,:) = tab_match_optcheby4.t_it(:);
%% V-MATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 8/log_optcheby1_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby1 = maketable(logfiles);

TSOLVE3(3,:) = tab_match_optcheby1.t_solve(:);
TITER3(3,:) = tab_match_optcheby1.t_it(:);
%% V-MATCH-4L1-JAC-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 8/log_l1jacobi_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_l1jac = maketable(logfiles);

TSOLVE3(4,:) = tab_match_l1jac.t_solve(:);
TITER3(4,:) = tab_match_l1jac.t_it(:);
%% Plot analysis
figure(1)
subplot(3,2,3)
loglog(NGPUS,TSOLVE3(1,:), ...
    NGPUS,TSOLVE3(2,:), ...
    NGPUS,TSOLVE3(3,:), ...
    NGPUS,TSOLVE3(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
figure(2)
subplot(3,2,3)
loglog(NGPUS,TITER3(1,:), ...
    NGPUS,TITER3(2,:), ...
    NGPUS,TITER3(3,:), ...
    NGPUS,TITER3(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on

%% V-MATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 10/log_cheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_cheby4 = maketable(logfiles);

TSOLVE4(1,:) = tab_match_cheby4.t_solve(:);
TITER4(1,:) = tab_match_cheby4.t_it(:);
NGPUS = tab_match_cheby4.np(:);
%% V-MATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 10/log_optcheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby4 = maketable(logfiles);

TSOLVE4(2,:) = tab_match_optcheby4.t_solve(:);
TITER4(2,:) = tab_match_optcheby4.t_it(:);
%% V-MATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 10/log_optcheby1_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby1 = maketable(logfiles);

TSOLVE4(3,:) = tab_match_optcheby1.t_solve(:);
TITER4(3,:) = tab_match_optcheby1.t_it(:);
%% V-MATCH-4L1-JAC-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 10/log_l1jacobi_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_l1jac = maketable(logfiles);

TSOLVE4(4,:) = tab_match_l1jac.t_solve(:);
TITER4(4,:) = tab_match_l1jac.t_it(:);
%% Plot analysis
figure(1)
subplot(3,2,4)
loglog(NGPUS,TSOLVE4(1,:), ...
    NGPUS,TSOLVE4(2,:), ...
    NGPUS,TSOLVE4(3,:), ...
    NGPUS,TSOLVE4(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
figure(2)
subplot(3,2,4)
loglog(NGPUS,TITER4(1,:), ...
    NGPUS,TITER4(2,:), ...
    NGPUS,TITER4(3,:), ...
    NGPUS,TITER4(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on

%% V-MATCH-4CHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 12/log_cheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_cheby4 = maketable(logfiles);

TSOLVE5(1,:) = tab_match_cheby4.t_solve(:);
TITER5(1,:) = tab_match_cheby4.t_it(:);
NGPUS = tab_match_cheby4.np(:);
%% V-MATCH-4OPTCHEBY4-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 12/log_optcheby4_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby4 = maketable(logfiles);

TSOLVE5(2,:) = tab_match_optcheby4.t_solve(:);
TITER5(2,:) = tab_match_optcheby4.t_it(:);
%% V-MATCH-4OPTCHEBY1-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 12/log_optcheby1_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_optcheby1 = maketable(logfiles);

TSOLVE5(3,:) = tab_match_optcheby1.t_solve(:);
TITER5(3,:) = tab_match_optcheby1.t_it(:);
%% V-MATCH-4L1-JAC-30l1JAC CPU/GPU
[~,logfiles] = system('ls -v1 12/log_l1jacobi_match*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match_l1jac = maketable(logfiles);

TSOLVE5(4,:) = tab_match_l1jac.t_solve(:);
TITER5(4,:) = tab_match_l1jac.t_it(:);
%% Plot analysis
figure(1)
subplot(3,2,5)
loglog(NGPUS,TSOLVE5(1,:), ...
    NGPUS,TSOLVE5(2,:), ...
    NGPUS,TSOLVE5(3,:), ...
    NGPUS,TSOLVE5(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
figure(2)
subplot(3,2,5)
loglog(NGPUS,TITER5(1,:), ...
    NGPUS,TITER5(2,:), ...
    NGPUS,TITER5(3,:), ...
    NGPUS,TITER5(4,:),'LineWidth',2)
xticks(NGPUS)
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on

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