% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturb the fixed points in the deterministic model to determine stability
% Case 2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% nexttile
% plot(Classes_ODE_S.t,Classes_ODE_S.S,'linewidth',2)
% hold on
% plot(Classes_ODE_Is.t,Classes_ODE_Is.S,'linewidth',2)
% hold on
% plot(Classes_ODE_Ia.t,Classes_ODE_Ia.S,'linewidth',2,'linestyle',':')
% hold on
% plot(Classes_ODE_Vo.t,Classes_ODE_Vo.S,'linewidth',2,'linestyle','--')
% xlabel('Time(days)')
% ylabel('Susceptibles')
% legend('S direction','Is direction','Ia direction')
% 
% nexttile
% plot(Classes_ODE_S.t,Classes_ODE_S.Is,'linewidth',2)
% hold on
% plot(Classes_ODE_Is.t,Classes_ODE_Is.Is,'linewidth',2)
% hold on
% plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Is,'linewidth',2,'linestyle',':')
% hold on
% plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Is,'linewidth',2,'linestyle','--')
% xlabel('Time(days)')
% ylabel('Symptomatic infections')
% legend('S direction','Is direction','Ia direction')
% 
% 
% nexttile
% plot(Classes_ODE_S.t,Classes_ODE_S.Ia,'linewidth',2)
% hold on
% plot(Classes_ODE_Is.t,Classes_ODE_Is.Ia,'linewidth',2)
% hold on
% plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Ia,'linewidth',2)
% hold on
% plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Ia,'linewidth',2)
% xlabel('Time(days)')
% ylabel('Asymptomatic infections')
% legend('S direction','Is direction','Ia direction')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Vipv,'linewidth',2)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Vipv,'linewidth',2,'linestyle','--')
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Vipv,'linewidth',2)
xlabel('Time(days)')
ylabel('IPV vaccinations')
legend('S direction','Is direction','Ia direction', 'Vopv direction')

% nexttile
% plot(Classes_ODE_S.t,Classes_ODE_S.Vopv,'linewidth',2)
% hold on
% plot(Classes_ODE_Is.t,Classes_ODE_Is.Vopv,'linewidth',2)
% hold on
% plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Vopv,'linewidth',2)
% hold on
% plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Vopv,'linewidth',2)
% xlabel('Time(days)')
% ylabel('OPV infections')
% legend('S direction','Is direction','Ia direction', 'Vopv direction')

nexttile
plot(Classes_ODE_S.t,Classes_ODE_S.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Is.t,Classes_ODE_Is.Stil,'linewidth',2)
hold on
plot(Classes_ODE_Ia.t,Classes_ODE_Ia.Stil,'linewidth',2,'linestyle','--')
hold on
plot(Classes_ODE_Vo.t,Classes_ODE_Vo.Stil,'linewidth',2)
xlabel('Time(days)')
ylabel('Partially immune')
legend('S direction','Is direction','Ia direction', 'Vopv direction')

















% %k determines which fixed point in this class we are looking at (ie where
% %we are on the line of fixed points).
% %We perturb our initial conditions by m.
% k=10000;
% m=1000;
% ICs = struct('S',0,'Is',0,'Ia',m,'Vipv',k,'Vopv',0,'Stil',para.N-k-m);
% [Classes_ODE] = ODE_polio_model(para,ICs,maxtime);
% 
% %Plot a phase portrait of Vipv, Stil and the variable in which we are
% %perturbing the initial condition along with the line of fixed points in
% %the (Vipv, Stil) plane and plot the starting position to make the
% %direction of motion clear.
% figure(3)
% clf
% plot3(Classes_ODE.Ia,Classes_ODE.Vipv,Classes_ODE.Stil)
% xlabel('Test variable')
% ylabel('Vipv')
% zlabel('Stil')
% hold on
% plot3(ICs.Ia,ICs.Vipv,ICs.Stil, '.')
% hold on
% %Plot the line of fixed points
% z=[1000:100:134000];
% plot3(zeros(length(z)),z,para.N-z);
% 
% % This analysis shows that given many different kinds of perturbation in all possible
% % directions, the trajectory moves from its initial position to some point
% % on the line of fixed points, indicating that this line of fixed points is
% % stable. Stability of individual points on this line is harder to
% % determine.
