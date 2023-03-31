%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot histograms side by side comparing control vs no control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note - this code only runs the control simulation - needs 'Polio_plots' to
%run first to compare against no control

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',14)

%Run model with contrl - change parameter values and the model algorithm
%being used to move between control schemes - currently set up to compare
%increased vaccination rate against no control
para_Control = struct('beta',1/7,'nu',0.0605,'gamma',1/42,'a',0.0005,'k',0.01,'q',1/21,'delta',1/60,'N',8799723);
m=0.05*para_Control.N;
ICs = struct('S',(1-m/para.N)*(1-0.904)*para.N,'Is',0,'Ia',m,'Vipv',(1-m/para_Control.N)*(1-0.619)*0.904*para.N,'Vopv',(1-m/para.N)*(1-0.105)*0.619*0.904*para.N,'Stil',(1-m/para.N)*0.105*0.619*0.904*para.N);

[Classes_Control] = Tauleap_polio_model(para_Control,ICs,maxtime,timestep);
NRuns=2000; 
IsCountControlMat = zeros(NRuns,length(Classes_Control.IsCount));

tic
for i=1:NRuns
    [Classes_Control] = Tauleap_polio_model(para_Control,ICs,maxtime,timestep);
    IsCountControlMat(i,:) = Classes_Control.IsCount;
end
toc

IsCount_Control_median = quantile(IsCountControlMat,0.5,1);

% %Plot two histograms on separate axes - unused in dissertation
% figure(1)
% clf
% tiledlayout(1,2)
% 
% nexttile
% %Plot the distribution of the final number paralysed individuals - no
% %control.
% histogram(IsCountMat(:,end),12,'Numbins',18,'facecolor','b','facealpha',0.2,'edgecolor','b','edgealpha',0.4,'linewidth',1)
% xlabel('Number of cases of paralysis')
% ylabel('Frequency')
% xline(IsCount_median(end),'b','linewidth',3)
% %xline(Classes_ODE.IsCount(end),'r','linewidth',3)
% legend('','Median','Deterministic')
% title('No Control')
% ax = gca; 
% ax.FontSize = f;
% 
% nexttile
% %Plot the distribution of the final number paralysed individuals - with
% %control
% histogram(IsCountControlMat(:,end),18,'facecolor','b','facealpha',0.2,'edgecolor','b','edgealpha',0.4,'linewidth',1)
% xlabel('Number of cases of paralysis')
% ylabel('Frequency')
% xline(IsCount_Control_median(end),'b','linewidth',3)
% %xline(Classes_Control_ODE.IsCount(end),'r','linewidth',3)
% legend('','Median','Deterministic')
% title('Increased vaccination rate of \nu = 0.0605')
% ax = gca; 
% ax.FontSize = f;

%Plot two histograms on same axes of different colours for different
%strategies
figure(2)
clf
hold on
histogram(IsCountMat(:,end),12,'Numbins',18,'facecolor','b','facealpha',0.2,'edgecolor','b','edgealpha',0.4,'linewidth',1)
histogram(IsCountControlMat(:,end),12,'Numbins',20,'facecolor',[0.8510,0.3255,0.0980],'facealpha',0.2,'edgecolor',[0.8510,0.3255,0.0980],'edgealpha',0.4,'linewidth',1)
xlabel('Number of cases of paralysis')
ylabel('Frequency')
xline(IsCount_median(end),'b','linewidth',3)
xline(IsCount_Control_median(end),'color',[0.8510,0.3255,0.0980],'linewidth',3)
xline(Classes_ODE.IsCount(end),'r','linewidth',3)
legend('','','No Control','Increased vaccination rate of \nu = 0.605','location','northeast')
%title('Quarantine of $I_a$','Interpreter','latex')
ax = gca; 
