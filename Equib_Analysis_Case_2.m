% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturb the fixed points in the deterministic model to determine stability
% Case 2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontsize',16)

%Choose parameters so the vaccine free equilibrium doesn't exist
para = struct('beta',1/7,'nu',0.0055,'gamma',1/42,'a',0.005,'k',0.01,'delta',1/60,'N',8799723);

%Identical analysis to before (but with the re-defined parameter values
%above)
maxtime = 1000;
%k selects which point along the line of fixed points we look at. Expect
%threshold at k = 1,392,649
k=1000000; 
%m is again the perturbation.
m=10;
ICsS = struct('S',m,'Is',0,'Ia',0,'Vipv',k-m/2,'Vopv',0,'Stil',para.N-k-m/2);
ICsIs = struct('S',0,'Is',m,'Ia',0,'Vipv',k-m/2,'Vopv',0,'Stil',para.N-k-m/2);
ICsIa = struct('S',0,'Is',0,'Ia',m,'Vipv',k-m/2,'Vopv',0,'Stil',para.N-k-m/2);
ICsVo = struct('S',0,'Is',0,'Ia',0,'Vipv',k-m/2,'Vopv',m,'Stil',para.N-k-m/2);
[Classes_ODE_S] = ODE_polio_model(para,ICsS,maxtime);
[Classes_ODE_Is] = ODE_polio_model(para,ICsIs,maxtime);
[Classes_ODE_Ia] = ODE_polio_model(para,ICsIa,maxtime);
[Classes_ODE_Vo] = ODE_polio_model(para,ICsVo,maxtime);


%Plot how Vipv and Stil change over time starting from perturbation -
%other classes investigated but unused in dissertation
figure(1)
clf
tiledlayout(1,2)

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Vipv,'linewidth',2,'linestyle','--')
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Vipv,'linewidth',2)
xlabel('Time (days)')
ylabel('$V_{IPV}$','interpreter','latex')
legend('$S$ direction','$I_s$ direction','$I_a$ direction', '$V_{OPV}$ direction','interpreter','latex','location','east')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Stil,'linewidth',2,'linestyle','--')
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Stil,'linewidth',2)
xlabel('Time (days)')
ylabel('$\widetilde{S}$','interpreter','latex')
legend('$S$ direction','$I_s$ direction','$I_a$ direction', '$V_{OPV}$ direction','interpreter','latex','location','east')

