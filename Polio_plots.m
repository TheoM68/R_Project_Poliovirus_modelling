% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run NRuns simulations - plot median of every class and plot shaded
% polygons for Is as well as the final distributions of 
% paralysed individuals (IsCount).
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',14)

%Set length of simulation, timestep and parameter values
maxtime = 365;
timestep = 1;
para = struct('beta',1/7,'nu',0.0055,'gamma',1/42,'a',0.0005,'k',0.01,'delta',1/60,'N',8799723);

%Define initial condition with m initial asymptomatic infections
m=0.05*para.N;
ICs = struct('S',(1-m/para.N)*(1-0.904)*para.N,'Is',0,'Ia',m,'Vipv',(1-m/para.N)*(1-0.619)*0.904*para.N,'Vopv',(1-m/para.N)*(1-0.105)*0.619*0.904*para.N,'Stil',(1-m/para.N)*0.105*0.619*0.904*para.N);

%Run initial simulations of stochastic and deterministic model.
[Classes] = Tauleap_polio_model(para,ICs,maxtime,timestep);
[Classes_ODE] = ODE_polio_model(para,ICs,maxtime);

%Define how many runs of our simulation we will do as well as matrices to
%store each class
NRuns=2000;                  
SMat = zeros(NRuns,length(Classes.S));
IsMat = zeros(NRuns,length(Classes.S));
IaMat = zeros(NRuns,length(Classes.S));
VopvMat = zeros(NRuns,length(Classes.S));
VipvMat = zeros(NRuns,length(Classes.S));
StilMat=zeros(NRuns,length(Classes.S));
IsCountMat = zeros(NRuns,length(Classes.S));

%Run the simulation NRuns times and store each class in each iteration as a
%row of its corresponding matrix.
tic
for i=1:NRuns
    [Classes] = Tauleap_polio_model(para,ICs,maxtime,timestep);
    SMat(i,:) = Classes.S;
    IsMat(i,:) = Classes.Is;
    IaMat(i,:) = Classes.Ia;
    VopvMat(i,:) = Classes.Vopv;
    VipvMat(i,:) = Classes.Vipv;
    StilMat(i,:) = Classes.Stil;
    IsCountMat(i,:) = Classes.IsCount;
end
toc

%Compute median value of each column (each day).
S_median = quantile(SMat,0.5,1);
Is_median = quantile(IsMat,0.5,1);
Ia_median = quantile(IaMat,0.5,1);
Vopv_median = quantile(VopvMat,0.5,1);
Vipv_median = quantile(VipvMat,0.5,1);
Stil_median = quantile(StilMat,0.5,1);
IsCount_median = quantile(IsCountMat,0.5,1);

%Calculates the 2.5% and 97.5% quantile and the median for each column of
%SMat - confidence intervals produced are only visible for Is so all of
%these computations excpet for Is are unused in report.
S_lower = quantile(SMat,0.025,1);
S_upper = quantile(SMat,0.975,1);
%Calculates the 2.5% and 97.5% quantile for each column of
%IsMat.
Is_lower = quantile(IsMat,0.025,1);
Is_upper = quantile(IsMat,0.975,1);
%Calculates the 2.5% and 97.5% quantile and the median for each column of
%IaMat.
Ia_lower = quantile(IaMat,0.025,1);
Ia_upper = quantile(IaMat,0.975,1);
%Calculates the 2.5% and 97.5% quantile and the median for each column of
%VopvMat.
Vopv_lower = quantile(VopvMat,0.025,1);
Vopv_upper = quantile(VopvMat,0.975,1);
%Calculates the 2.5% and 97.5% quantile and the median for each column of
%VipvMat.
Vipv_lower = quantile(VipvMat,0.025,1);
Vipv_upper = quantile(VipvMat,0.975,1);
%Calculates the 2.5% and 97.5% quantile and the median for each column of
%StilMat.
Stil_lower = quantile(StilMat,0.025,1);
Stil_upper = quantile(StilMat,0.975,1);
%Calculates the 2.5% and 97.5% quantile for each column of
%IsCountMat.
IsCount_lower = quantile(IsCountMat,0.025,1);
IsCount_upper = quantile(IsCountMat,0.975,1);

%Define edges of shaded polygon for 95% confidence - all but Is unused in report
t_poly = [Classes.t fliplr(Classes.t)];
S_poly = [S_upper fliplr(S_lower)];
Is_poly = [Is_upper fliplr(Is_lower)];
Ia_poly = [Ia_upper fliplr(Ia_lower)];
Vopv_poly = [Vopv_upper fliplr(Vopv_lower)];
Vipv_poly = [Vipv_upper fliplr(Vipv_lower)];
Stil_poly = [Stil_upper fliplr(Stil_lower)];
IsCount_poly = [IsCount_upper fliplr(IsCount_lower)];

%Total force of infection term to compute final size
for i = 1:NRuns
    FOI(i,:) = para.beta/para.N*(IaMat(i,:) + IsMat(i,:)).*SMat(i,:) + ...
    para.beta/para.N*(IaMat(i,:) + IsMat(i,:)).*VipvMat(i,:) + ...
    para.k*para.beta/para.N*(IaMat(i,:) + IsMat(i,:)).*VopvMat(i,:) + ...
    para.k*para.beta/para.N*(IaMat(i,:) + IsMat(i,:)).*StilMat(i,:) + para.delta*VopvMat(i,:);
end
%Sum all force of infection terms to get final size - shows no failed outbreaks
FinalSize = sum(FOI,2);

%Colour mapping slightly darker than default colours to show determinstic
%solution on top of median stochastic dynamics
CMapDark = [0,0.3176,0.5294;...
    0.7098,0.2471,0.0510;...
    0.7804,0.5294,0.0627;...
    0,0,0;...
    0.3216,0.4902,0.0980;...
    0.1647,0.5255,0.6784];

%Tiled plot to show dyanimcs of whole system, dynamics of Is and histogram
%of paralytic cases
figure(1)
clf
t = tiledlayout(1,3,'TileSpacing','Compact');

nexttile
%Plot the median of each class on the same axes
plot(Classes.t,S_median,Classes.t,Is_median,Classes.t,Ia_median,Classes.t,Vopv_median,Classes.t,Vipv_median,Classes.t,Stil_median,'linewidth',2)
% %Plot Determinstic solutions with darker colours - only for inital plot
% hold on
% plot(Classes_ODE.t,Classes_ODE.S,'--','color',CMapDark(1,:))
% plot(Classes_ODE.t,Classes_ODE.Is,'--','color',CMapDark(2,:))
% plot(Classes_ODE.t,Classes_ODE.Ia,'--','color',CMapDark(3,:))
% plot(Classes_ODE.t,Classes_ODE.Vopv,'--','color',CMapDark(4,:))
% plot(Classes_ODE.t,Classes_ODE.Vipv,'--','color',CMapDark(5,:))
% plot(Classes_ODE.t,Classes_ODE.Stil,'--','color',CMapDark(6,:))
xlabel('Time (days)')
ylabel('Number of individuals')
legend('$S$','$I_s$','$I_a$','$V_{OPV}$','$V_{IPV}$','$\widetilde{S}$','location','east','interpreter','latex')
axis([0 365 0 inf])

nexttile
%Plot the shaded polygon showing 95% confidence region for Is
fill(t_poly,Is_poly,'b','facealpha',0.2,'edgecolor','none')
xlabel('Time (days)')
ylabel('$I_s$','interpreter','latex')
hold on
%Add median line on top
plot(Classes.t,Is_median,'b','Linewidth',2)
hold on
plot(Classes_ODE.t,Classes_ODE.Is,'r','Linewidth',1.5)
legend('2.5th and 97.5th percentiles','Median','Deterministic')
axis([0 365 0 inf])

nexttile
%Plot the distribution of the final number paralysed individuals
histogram(IsCountMat(:,end),12,'Numbins',18,'facecolor','b','facealpha',0.2,'edgecolor','b','edgealpha',0.4,'linewidth',1)
xlabel('Number of cases of paralysis')
ylabel('Frequency')
xline(IsCount_median(end),'b','linewidth',3)
xline(Classes_ODE.IsCount(end),'r','linewidth',3)
legend('','Median','Deterministic')

t.Padding = 'none';

set(gcf,'windowstyle','normal')
set(gcf,'position',[343,439,1377,458])




% %Tiled plot for deterministic output of whole system and stochastic - with
% %confidence intervals - they are impossible to tell apart so unused in
% %dissertation
% figure(2)
% clf
% hold on
% plot(Classes.t,S_median,Classes.t,Is_median,Classes.t,Ia_median,Classes.t,Vopv_median,Classes.t,Vipv_median,Classes.t,Stil_median,'linewidth',2)
% fill(t_poly,S_poly,'b','facealpha',0.2,'edgecolor','none')
% fill(t_poly,Is_poly,'b','facealpha',0.2,'edgecolor','none')
% fill(t_poly,Ia_poly,'b','facealpha',0.2,'edgecolor','none')
% fill(t_poly,Vopv_poly,'b','facealpha',0.2,'edgecolor','none')
% fill(t_poly,Vipv_poly,'b','facealpha',0.2,'edgecolor','none')
% fill(t_poly,Stil_poly,'b','facealpha',0.2,'edgecolor','none')
% plot(Classes_ODE.t,Classes_ODE.S,Classes_ODE.t,Classes_ODE.Is,Classes_ODE.t,Classes_ODE.Ia,Classes_ODE.t,Classes_ODE.Vopv,Classes_ODE.t,Classes_ODE.Vipv,Classes_ODE.t,Classes_ODE.Stil,'linewidth',2)
% xlabel('Time(days)')
% ylabel('Number of individuals')
% legend('S','Is','Ia','Vopv','Vipv','Stil','location','east')
% axis([0 365 0 inf])
% ax = gca; 
% ax.FontSize = f;


% 
% %Plot the shaded polygon showing 2.5th and 97.5th percentiles for IsCoun - unused
% figure(3)
% clf
% fill(t_poly,IsCount_poly,'b','facealpha',0.2,'edgecolor','none')
% xlabel('Time(days)')
% ylabel('Individuals with paralysis')
% hold on
% %Add median line on top
% plot(Classes.t,IsCount_median,'b','Linewidth',2)
% hold on
% plot(Classes_ODE.t,Classes_ODE.IsCount,'r','Linewidth',1.5)









