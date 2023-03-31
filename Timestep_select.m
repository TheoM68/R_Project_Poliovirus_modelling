% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot histograms of how many times our tau-leaping algorithm produces
% negative populations for 3 different time steps.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',14) 

%Set up inputs for model as described in model inputs in model report
maxtime = 365;
para = struct('beta',1/7,'nu',0.0055,'gamma',1/42,'a',0.005,'k',0.01,'delta',1/60,'N',8799723);

m=0.05*para.N;
ICs = struct('S',(1-m/para.N)*(1-para.nu)*para.N,'Is',0,'Ia',m,'Vipv',(1-m/para.N)*(1-0.619)*para.nu*para.N,'Vopv',(1-m/para.N)*(1-0.105)*0.619*para.nu*para.N,'Stil',(1-m/para.N)*0.105*0.619*para.nu*para.N);

%Set number of runs of our model
NRuns = 2000;

%Set first timestep and matrix for storing 'negative count'
timestep = 50;
A = zeros([NRuns,1]);
%Run our model NRuns times with first timestep and store final 'negative count'
for i=1:NRuns
    [Classes] = Tauleap_polio_model(para,ICs,maxtime,timestep);
    A(i)=Classes.NegativeCount(end);
end
%Calculate the median final 'negative count' across all simulations
A_median = quantile(A,0.5,1);

%Rerun for timestep of 10
timestep = 10;
B = zeros([NRuns,1]);

for i=1:NRuns
    [Classes] = Tauleap_polio_model(para,ICs,maxtime,timestep);
    B(i)=Classes.NegativeCount(end);
end

B_median = quantile(B,0.5,1);

%Rerun for timestep of 1
timestep = 1;
C = zeros([NRuns,1]);

for i=1:NRuns
    [Classes] = Tauleap_polio_model(para,ICs,maxtime,timestep);
    C(i)=Classes.NegativeCount(end);
end

C_median = quantile(C,0.5,1);

%Set color of historgrams
col = [0.4 0 1];

%Plot a histogram for each timestep
figure(1)
clf
t = tiledlayout(1,3,'Tilespacing','Compact');

nexttile
histogram(A,'facecolor',col,'facealpha',0.2,'edgecolor',col,'edgealpha',0.4,'linewidth',1)
xA=xline(A_median,'linewidth',3);
xlabel('Negative count')
set(xA,'color',col)
ylabel('Frequency')
title('\tau = 50')
axis([-1 8 0 1087])

nexttile
histogram(B,'facecolor',col,'facealpha',0.2,'edgecolor',col,'edgealpha',0.4,'linewidth',1)
xB=xline(B_median,'linewidth',3);
set(xB,'color',col)
xlabel('Negative count')
ylabel('Frequency')
title('\tau = 10')
axis([-1 8 0 1981])

nexttile
histogram(C,'facecolor',col,'facealpha',0.2,'edgecolor',col,'edgealpha',0.4,'linewidth',1)
xC=xline(C_median,'linewidth',3);
set(xC,'color',col)
xlabel('Negative count')
ylabel('Frequency')
title('\tau = 1')
axis([-1 8 0 2001])
legend(['Distribution of' newline 'negative count.'],'Median','location','east')

t.Padding = 'none';

set(gcf,'windowstyle','normal')
set(gcf,'position',[343,439,1377,458])

