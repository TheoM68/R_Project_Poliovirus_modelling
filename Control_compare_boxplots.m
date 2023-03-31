%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare distribution of paralysis cases of all control schemes using 
% box plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontsize',14)

%Create a matrix with columns as the final distributions of
%paralytic cases across 200 simulations for each control scheme.

maxtime = 365;
timestep = 1;

%Define initial condition with m initial asymptomatic infections
m=0.05*para1.N;
ICs = struct('S',(1-m/para1.N)*(1-0.904)*para1.N,'Is',0,'Ia',m,'Vipv',(1-m/para1.N)*(1-0.619)*0.904*para1.N,'Vopv',(1-m/para1.N)*(1-0.105)*0.619*0.904*para1.N,'Stil',(1-m/para1.N)*0.105*0.619*0.904*para1.N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%No control - base case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para1 = struct('beta',1/7,'nu',0.0055,'gamma',1/42,'a',0.0005,'k',0.01,'delta',1/60,'N',8799723);

%Base model with no control to compare other models against.
[Classes_base] = Tauleap_polio_model(para1,ICs,maxtime,timestep);

%Matrix of final paralysis cases
NRuns=2000;                  
IsCountBaseMat = zeros(NRuns,length(Classes_base.Is));

tic
for i=1:NRuns
    [Classes_base] = Tauleap_polio_model(para1,ICs,maxtime,timestep);
    IsCountBaseMat(i,:) = Classes_base.IsCount;
end
toc

IsCountBaseFinal = IsCountBaseMat(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control implementing higher vaccination rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para2 = struct('beta',1/7,'nu',0.0605,'gamma',1/42,'a',0.0005,'k',0.01,'delta',1/60,'N',8799723);

[Classes_vaccine] = Tauleap_polio_model(para2,ICs,maxtime,timestep);

%Matrix of final paralysis cases
NRuns=2000;                  
IsCountVaccineMat = zeros(NRuns,length(Classes_vaccine.Is));

tic
for i=1:NRuns
    [Classes_vaccine] = Tauleap_polio_model(para2,ICs,maxtime,timestep);
    IsCountVaccineMat(i,:) = Classes_vaccine.IsCount;
end
toc

IsCountVaccineFinal = IsCountVaccineMat(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control implementing a reduced contact rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para3 = struct('beta',1/14,'nu',0.0055,'gamma',1/42,'a',0.0005,'k',0.01,'delta',1/60,'N',8799723);

[Classes_contact] = Tauleap_polio_model(para3,ICs,maxtime,timestep);

%Matrix of final paralysis cases
NRuns=2000;                  
IsCountContactMat = zeros(NRuns,length(Classes_contact.Is));

tic
for i=1:NRuns
    [Classes_contact] = Tauleap_polio_model(para3,ICs,maxtime,timestep);
    IsCountContactMat(i,:) = Classes_contact.IsCount;
end
toc

IsCountContactFinal = IsCountContactMat(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control implementing quarantine of Ia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para4 = struct('beta',1/7,'nu',0.0055,'gamma',1/42,'a',0.0005,'k',0.01,'delta',1/60,'q',1/21,'N',8799723);

[Classes_quar_Ia] = Tauleap_polio_model_quarantine_Ia(para4,ICs,maxtime,timestep);

%Matrix of final paralysis cases
NRuns=2000;                  
IsCountQuarIaMat = zeros(NRuns,length(Classes_quar_Ia.Is));

tic
for i=1:NRuns
    [Classes_quar_Ia] = Tauleap_polio_model_quarantine_Ia(para4,ICs,maxtime,timestep);
    IsCountQuarIaMat(i,:) = Classes_quar_Ia.IsCount;
end
toc

IsCountQuarIaFinal = IsCountQuarIaMat(:,end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control implementing quarantine of Is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para5 = struct('beta',1/7,'nu',0.0055,'gamma',1/42,'a',0.0005,'k',0.01,'delta',1/60,'q',1/21,'N',8799723);

[Classes_quar_Is] = Tauleap_polio_model_quarantine_Is(para5,ICs,maxtime,timestep);

%Matrix of final paralysis cases
NRuns=2000;                  
IsCountQuarIsMat = zeros(NRuns,length(Classes_quar_Is.Is));

tic
for i=1:NRuns
    [Classes_quar_Is] = Tauleap_polio_model_quarantine_Is(para4,ICs,maxtime,timestep);
    IsCountQuarIsMat(i,:) = Classes_quar_Is.IsCount;
end
toc

IsCountQuarIsFinal = IsCountQuarIsMat(:,end);

%Matrix of cases of paralysis with each column a control scheme in order of
%effectiveness
IsCountCompare = [IsCountQuarIsFinal,IsCountBaseFinal,IsCountContactFinal,IsCountQuarIaFinal,IsCountVaccineFinal];

%Determine colour of box plots
col = [0 0.5216 0.9804];

%Plot box plot
figure(1)
clf
p1 = boxplot(IsCountCompare,'labels',{'Quarantine of $I_s$','No intervention','Reduced contact rate','Quarantine of $I_a$','Vaccine booster'},...
    'Colors',[0.2314,0.2824,1],'symbol','+');
ylabel("Number of paralytic cases")

%Recolour and resize all lines
for ih=1:6
    set(p1(ih,:),'LineWidth',1.6,'color',col);
end
%Recolour the median lines
medlines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(medlines, 'Color',[0.0549 0.4039 0.7098]);

%Change linewidth of outliers
set(p1(7,:),'LineWidth',1.2);
%Size and colour of outliers
outl = findobj(gcf,'tag','Outliers');
for j = 1:length(outl)
   outl(j).MarkerEdgeColor = [1 1 1];
   outl(j).MarkerSize = 5;
end

%Fill boxes 
box_h = findobj(p1,'Tag','Box');
for i = 1:length(box_h)
    hp(i) = patch([box_h(i).XData],[box_h(i).YData],col,'FaceAlpha',0.3,'LineStyle','none');
end

gc = gca;
%Ensure boxplot plotted on top of coloured patches
gc.Children = gc.Children([6 1 2 3 4 5]);
%Format x labels properly
gc.XAxis.TickLabelInterpreter = 'latex';

%Move upper whisker to 97.5th percentile
a=findobj(gcf,'Tag','Upper Adjacent Value');
a(1).YData = [quantile(IsCountVaccineFinal,0.975,1) quantile(IsCountVaccineFinal,0.975,1)];
a(2).YData = [quantile(IsCountQuarIaFinal,0.975,1) quantile(IsCountQuarIaFinal,0.975,1)];
a(3).YData = [quantile(IsCountContactFinal,0.975,1) quantile(IsCountContactFinal,0.975,1)];
a(4).YData = [quantile(IsCountBaseFinal,0.975,1) quantile(IsCountBaseFinal,0.975,1)];
a(5).YData = [quantile(IsCountQuarIsFinal,0.975,1) quantile(IsCountQuarIsFinal,0.975,1)];

b=findobj(gcf,'Tag','Upper Whisker');
b(1).YData(2) = quantile(IsCountVaccineFinal,0.975,1);
b(2).YData(2) = quantile(IsCountQuarIaFinal,0.975,1);
b(3).YData(2) = quantile(IsCountContactFinal,0.975,1);
b(4).YData(2) = quantile(IsCountBaseFinal,0.975,1);
b(5).YData(2) = quantile(IsCountQuarIsFinal,0.975,1);

%Move lower whisker to 2.5th percentile
c=findobj(gcf,'Tag','Lower Adjacent Value');
c(1).YData = [quantile(IsCountVaccineFinal,0.025,1) quantile(IsCountVaccineFinal,0.025,1)];
c(2).YData = [quantile(IsCountQuarIaFinal,0.025,1) quantile(IsCountQuarIaFinal,0.025,1)];
c(3).YData = [quantile(IsCountContactFinal,0.025,1) quantile(IsCountContactFinal,0.025,1)];
c(4).YData = [quantile(IsCountBaseFinal,0.025,1) quantile(IsCountBaseFinal,0.025,1)];
c(5).YData = [quantile(IsCountQuarIsFinal,0.025,1) quantile(IsCountQuarIsFinal,0.025,1)];

d=findobj(gcf,'Tag','Lower Whisker');
d(1).YData(1) = quantile(IsCountVaccineFinal,0.025,1);
d(2).YData(1) = quantile(IsCountQuarIaFinal,0.025,1);
d(3).YData(1) = quantile(IsCountContactFinal,0.025,1);
d(4).YData(1) = quantile(IsCountBaseFinal,0.025,1);
d(5).YData(1) = quantile(IsCountQuarIsFinal,0.025,1);
