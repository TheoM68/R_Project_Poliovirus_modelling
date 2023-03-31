% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot effective reproductive number
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set length of simulation, timestep and parameter values
maxtime = 365;
timestep = 1;
para = struct('beta',1/7,'nu',0.904,'gamma',1/42,'a',0.005,'k',0.01,'delta',1/60,'N',8799723);

%Define initial condition with m initial asymptomatic infections
m=0.05*para.N;
ICs = struct('S',(1-m/para.N)*(1-para.nu)*para.N,'Is',0,'Ia',m,'Vipv',(1-m/para.N)*(1-0.619)*para.nu*para.N,'Vopv',(1-m/para.N)*(1-0.105)*0.619*para.nu*para.N,'Stil',(1-m/para.N)*0.105*0.619*para.nu*para.N);

%Run model
[Classes] = ODE_polio_model(para,ICs,maxtime);

%Compute Re
Re = para.beta/(para.gamma*para.N)*(Classes.S + Classes.Vipv + para.k.*Classes.Vopv + para.k*Classes.Stil);

%Plot Re as a function of time
figure(1)
clf
plot(Classes.t,Re,'linewidth',2)
axis([0 365 0 inf])
xlabel('Time (days)')
ylabel('Effective reproductive number R_e(t)')
ax = gca; 
ax.FontSize = 12;