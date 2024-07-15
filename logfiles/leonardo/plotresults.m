clear all; close all
load data2.mat

%3d Laplacian matching-based aggregation

celle=tab_3dlapmatch;
n=size(celle,1);
it=[];
tsolve=[];
tit=[];
for i=1:n
    tab=celle{i,1}
   name=tab(1,1)
   % name=table2array(name)
   % methd=char(name)
   % methd=methd(11:end-8)
    
    %itname=matlab.lang.makeValidName(strcat('it',methd))
    %tsolvename=matlab.lang.makeValidName(strcat('tsolve',methd))
    %titname=matlab.lang.makeValidName(strcat('tit',methd))
    
    val=tab(:,2:end)
    valarray=table2array(val);

    np=valarray(:,1);
    [np,idx]=sort(np);
    
    %assignin('base',itname,valarray(idx,7))
    it=[it; valarray(idx,7)']
    %assignin('base',tsolvename,valarray(idx,6))
    tsolve=[tsolve; valarray(idx,6)']
    %assignin('base',titname,valarray(idx,8))
    tit=[tit; valarray(idx,8)']
    
end
numconf=length(np)

%% Plot analysis
%deg 4
idx4=[1,6,11,16]
figure(1)
subplot(3,2,1)
bar([1:numconf],tsolve(idx4,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,1)
bar([1:numconf],tit(idx4,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,1)
bar([1:numconf],it(idx4,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%
%% Plot analysis
%deg 6
idx6=[2,7,12,17]
figure(1)
subplot(3,2,2)
bar([1:numconf],tsolve(idx6,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,2)
bar([1:numconf],tit(idx6,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,2)
bar([1:numconf],it(idx6,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%
%% Plot analysis
%deg 8
idx8=[3,8,13,18]
figure(1)
subplot(3,2,3)
bar([1:numconf],tsolve(idx8,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,3)
bar([1:numconf],tit(idx8,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,3)
bar([1:numconf],it(idx8,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%% Plot analysis
%deg 10
idx10=[4,9,14,19]
figure(1)
subplot(3,2,4)
bar([1:numconf],tsolve(idx10,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,4)
bar([1:numconf],tit(idx10,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,4)
bar([1:numconf],it(idx10,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%% Plot analysis
%deg 12
idx12=[5,10,15,20]
figure(1)
subplot(3,2,5)
bar([1:numconf],tsolve(idx12,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,5)
bar([1:numconf],tit(idx12,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,5)
bar([1:numconf],it(idx12,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%
pause

keyboard()

print(figure(1), '-djpeg', '3dlapmatch-tsolve.jpg');
print(figure(2), '-djpeg', '3dlapmatch-tit.jpg');
print(figure(3), '-djpeg', '3dlapmatch-it.jpg');


close all
pause

%3d laplacian, VBM aggregation

celle=tab_3dlapsoc1;
n=size(celle,1);
it=[];
tsolve=[];
tit=[];
for i=1:n
    tab=celle{i,1}
    name=tab(1,1)
    %name=table2array(name)
    %methd=char(name)
    %methd=methd(11:end-8)
    
    %itname=matlab.lang.makeValidName(strcat('it',methd))
    %tsolvename=matlab.lang.makeValidName(strcat('tsolve',methd))
    %titname=matlab.lang.makeValidName(strcat('tit',methd))
    
    val=tab(:,2:end)
    valarray=table2array(val);

    np=valarray(:,1);
    [np,idx]=sort(np);
    
    %assignin('base',itname,valarray(idx,7))
    it=[it; valarray(idx,7)']
    %assignin('base',tsolvename,valarray(idx,6))
    tsolve=[tsolve; valarray(idx,6)']
    %assignin('base',titname,valarray(idx,8))
    tit=[tit; valarray(idx,8)']
    
end
numconf=length(np)

%% Plot analysis
%deg 4
idx4=[1,6,11,16]
figure(1)
subplot(3,2,1)
bar([1:numconf],tsolve(idx4,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,1)
bar([1:numconf],tit(idx4,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,1)
bar([1:numconf],it(idx4,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 4/4 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%
%% Plot analysis
%deg 6
idx6=[2,7,12,17]
figure(1)
subplot(3,2,2)
bar([1:numconf],tsolve(idx6,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,2)
bar([1:numconf],tit(idx6,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,2)
bar([1:numconf],it(idx6,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 6/6 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%
%% Plot analysis
%deg 8
idx8=[3,8,13,18]
figure(1)
subplot(3,2,3)
bar([1:numconf],tsolve(idx8,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,3)
bar([1:numconf],tit(idx8,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,3)
bar([1:numconf],it(idx8,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 8/8 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%% Plot analysis
%deg 10
idx10=[4,9,14,19]
figure(1)
subplot(3,2,4)
bar([1:numconf],tsolve(idx10,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,4)
bar([1:numconf],tit(idx10,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,4)
bar([1:numconf],it(idx10,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 10/10 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
%% Plot analysis
%deg 12
idx12=[5,10,15,20]
figure(1)
subplot(3,2,5)
bar([1:numconf],tsolve(idx12,:)')
set(gca,'XTickLabel',{'64','128','256','512','1024'});
xlabel('Number of GPUs')
ylabel('Solve time (s)')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(2)
subplot(3,2,5)
bar([1:numconf],tit(idx12,:)')
set(gca,'XTickLabel',{'64','256','128','512','1024'});
xlabel('Number of GPUs')
ylabel('Time per iteration (s)')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')
figure(3)
subplot(3,2,5)
bar([1:numconf],it(idx12,:)')
set(gca,'XTickLabel',{'64','128','256','512', '1024'});
xlabel('Number of GPUs')
ylabel('Iterations')
title('Polynomial of degree 12/12 L1-JACOBI sweeps')
axis tight
grid on
legend({'Chebyshev 4','Opt. Chebyshev 4', ...
    'Opt. Chebyshev 1','L1-Jacobi'},'Location','northwest')

pause

print(figure(1), '-djpeg', '3dlapsoc1-tsolve.jpg');
print(figure(2), '-djpeg', '3dlapsoc1-tit.jpg');
print(figure(3), '-djpeg', '3dlapsoc1-it.jpg');