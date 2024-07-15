%% Analysis of the results
% The script uses matlab2tikz to produce figures:
% https://github.com/matlab2tikz/matlab2tikz
% and some bash-command to parse the script, thus it can be run only on a
% Linux machine.

clear; clc; close all;
%addpath("../../matlab2tikz/src/"); % or where it is

%% SOC1
[~,logfiles] = system('ls -v1 soc1/log_soc1_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_soc1 = maketable(logfiles);

%% MATCH
[~,logfiles] = system('ls -v1 match/log_match_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_match = maketable(logfiles);

%% 3DLAP - SOC1
[~,logfiles] = system('ls -v1 3dlap/soc1/log_poly_soc1_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_3dlapsoc1 = maketable(logfiles);

%% 3DLAP - MATCH
[~,logfiles] = system('ls -v1 3dlap/match/log_poly_match_l1jac*');
logfiles = split(logfiles);
logfiles(end) = [];
tab_3dlapmatch = maketable(logfiles);

%% Save data to file
save('data2.mat',"tab_3dlapmatch","tab_3dlapsoc1","tab_match","tab_soc1");

%% Function to analyze logs
function tab = maketable(logfiles)
%MAKETABLE produces tables with the result of a given run
n = length(logfiles);

for i = 1:n
    file = logfiles{i};

    [~,out] = system(sprintf('grep "Building new smoother: " %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j = 1:smtnumber
        try
            smother{i,j} = sscanf(outsmoothers{j},'Building new smoother: %s');
        catch
            smother{i,j} = '';
        end
    end

    [~,out] = system(sprintf('grep "Iterations to convergence          :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j = 1:smtnumber
        try
            iter(i,j) = sscanf(outsmoothers{j},'Iterations to convergence          : %d');
        catch
            iter(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Time to build hierarchy            :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j = 1:smtnumber
        try
            thierarchy(i,j) = sscanf(outsmoothers{j},'Time to build hierarchy            : %e');
        catch
            thierarchy(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Time to build smoothers            :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j = 1:smtnumber
        try
            tsmoothers(i,j) = sscanf(outsmoothers{j},'Time to build smoothers            : %e');
        catch
            tsmoothers(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Total time for preconditioner      :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j = 1:smtnumber
        try
            tbuild(i,j) = sscanf(outsmoothers{j},'Total time for preconditioner      : %e');
        catch
            tbuild(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Time to solve system               :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            tsolve(i,j) = sscanf(outsmoothers{j},'Time to solve system               : %e');
        catch
            tsolve(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Time per iteration                 :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            titer(i,j) = sscanf(outsmoothers{j},'Time per iteration                 : %e');
        catch
            titer(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Total time                         :" %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            tottime(i,j) = sscanf(outsmoothers{j},'Total time                         : %e');
        catch
            tottime(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Linear system size                 : " %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            sysize(i,j) = sscanf(outsmoothers{j},'Linear system size                 : %e');
        catch
            sysize(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "   Number of levels   :            " %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            levels(i,j) = sscanf(outsmoothers{j},'   Number of levels   : %d');
        catch
            levels(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "   Operator complexity:    " %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            opc(i,j) = sscanf(outsmoothers{j},'   Operator complexity: %f');
        catch
            opc(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "   Average coarsening :    " %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            avgcoars(i,j) = sscanf(outsmoothers{j},'   Average coarsening : %f');
        catch
            avgcoars(i,j) = NaN;
        end
    end

    [~,out] = system(sprintf('grep "Computed solution on " %s',file));
    outsmoothers = splitlines(out);
    outsmoothers(end) = [];
    smtnumber = length(outsmoothers);
    for j=1:smtnumber
        try
            np(i,j) = sscanf(outsmoothers{j},'Computed solution on %d');
        catch
            np(i,j) = NaN;
        end
    end
end

tab = cell(smtnumber,1);
for j = 1:smtnumber
    tab{j} = table(smother(:,j),np(:,j),sysize(:,j),thierarchy(:,j),tsmoothers(:,j), ...
    tbuild(:,j),tsolve(:,j),iter(:,j),titer(:,j), ...
    opc(:,j),avgcoars(:,j),levels(:,j),'VariableNames',{'prec','np','n','t_hier','t_smooth', ...
    't_build','t_solve','it','t_it','opc','avg_coars','lev'});
end

end