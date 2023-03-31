% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturb the fixed points in the deterministic model to determine stability
% Case 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set maxtime - greater than a year since we need longer to see the system
%equilibriate
maxtime = 1000;

%Investigate fixed point (1) (non-zero asymptomatic infections), perturb
%the Stil coordinate since this was the coordinate with neutral stability
%in the linear system.

%Parameters chosen based on data and beta changed to ensure fixed point 1
%exists
para = struct('beta',3,'nu',0.0055,'gamma',1/42,'a',0.005,'k',0.01,'delta',1/60,'N',8799723);

%Perturb by m and place m individuals into all other classes
m=1;
ICsS = struct('S',m,'Is',0,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',0,'Vopv',0,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
ICsIs = struct('S',0,'Is',m,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',0,'Vopv',0,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
ICsVi = struct('S',0,'Is',0,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',m,'Vopv',0,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
ICsVo = struct('S',0,'Is',0,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',0,'Vopv',m,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
[Classes_ODE_S] = ODE_polio_model(para,ICsS,maxtime);
[Classes_ODE_Is] = ODE_polio_model(para,ICsIs,maxtime);
[Classes_ODE_Vi] = ODE_polio_model(para,ICsVi,maxtime);
[Classes_ODE_Vo] = ODE_polio_model(para,ICsVo,maxtime);

% Plot asymptomatic infections and Stil against time for
% all 4 different perturbation directions and 2 different perturbations
figure(1)
clf
tiledlayout(2,2)

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Vi.t,Classes_ODE_Vi.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Ia,'linewidth',2)
xlabel('Time(days)')
ylabel('Asymptomatic infections')
legend('S direction','I_s direction','Vipv direction', 'Vopv direction')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Vi.t,Classes_ODE_Vi.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Stil,'linewidth',2)
xlabel('Time(days)')
ylabel('Partially immune')
legend('S direction','Is direction','Vipv direction', 'Vopv direction')

% Rerun simulation with much larger perturbation of 1000000.
m=1000000;
ICsS = struct('S',m,'Is',0,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',0,'Vopv',0,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
ICsIs = struct('S',0,'Is',m,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',0,'Vopv',0,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
ICsVi = struct('S',0,'Is',0,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',m,'Vopv',0,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
ICsVo = struct('S',0,'Is',0,'Ia',para.N - para.gamma*para.N/(para.k*para.beta)-m/2,'Vipv',0,'Vopv',m,'Stil',para.gamma*para.N/(para.k*para.beta)-m/2);
[Classes_ODE_S] = ODE_polio_model(para,ICsS,maxtime);
[Classes_ODE_Is] = ODE_polio_model(para,ICsIs,maxtime);
[Classes_ODE_Vi] = ODE_polio_model(para,ICsVi,maxtime);
[Classes_ODE_Vo] = ODE_polio_model(para,ICsVo,maxtime);

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Vi.t,Classes_ODE_Vi.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Ia,'linewidth',2)
xlabel('Time(days)')
ylabel('Asymptomatic infections')
legend('S direction','Is direction','Vipv direction', 'Vopv direction')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Vi.t,Classes_ODE_Vi.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Stil,'linewidth',2)
xlabel('Time(days)')
ylabel('Partially immune')
legend('S direction','Is direction','Vipv direction', 'Vopv direction')


%For these same parameter values where the above fixed point exists we
%investigate stability of the other class of fixed points (2) without
%infections 

%Extend simulation time to see the system equilibriate.
maxtime = 4000;
%k selects which point along the line of fixed points we look at.
k=100;
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

figure(2)
clf
tiledlayout(1,3)

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Ia,'linewidth',3)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Ia,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Ia,'linewidth',2)
xlabel('Time(days)')
ylabel('Asymptomatic infections')
legend('S direction','Is direction','Ia direction', 'Vopv direction','location','northwest')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Vipv,'linewidth',2)
xlabel('Time(days)')
ylabel('IPV vaccinations')
legend('S direction','Is direction','Ia direction', 'Vopv direction')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Stil,'linewidth',2)
xlabel('Time(days)')
ylabel('Partially immune')
legend('S direction','Is direction','Ia direction', 'Vopv direction')

%This whole class displays instability due to the threshold for Vipv for
%isntability being negative in this case hence it is always satisfied hence
%we always see instability.
