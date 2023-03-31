%Deterministic poliovirus model
function [Classes] = ODE_polio_model(para,ICs,maxtime)

    %Run ODE using ODE45

    opts = odeset('RelTol',1e-5);
    [t, pop] = ode45(@diff_polio_model, [0:1: maxtime], [ICs.S ICs.Is ICs.Ia ICs.Vipv ICs.Vopv ICs.Stil 0], opts, para);

    %Convert output to structure
    Classes = struct('S',pop(:,1),'Is',pop(:,2),'Ia',pop(:,3),'Vipv',pop(:,4),'Vopv',pop(:,5),'Stil',pop(:,6),'IsCount',pop(:,7),'t',t);


    %Diff equations
    function dPop = diff_polio_model(t,pop,para)

        S=pop(1);
        Is=pop(2);
        Ia=pop(3);
        Vipv=pop(4);
        Vopv=pop(5);
        Stil=pop(6);
        %Add a class which counts new symptomatic infections
        IsCount=pop(7);

        dS = -para.beta*S*(Ia + Is)/para.N - para.nu*S;
        dIs = para.a*para.beta*S*(Ia + Is)/para.N - para.gamma*Is;
        dIa = (para.beta*(Ia + Is)/para.N)*((1-para.a)*S + Vipv + para.k*Vopv + para.k*Stil) + para.delta*Vopv - para.gamma*Ia;
        dVipv = para.nu*S - para.beta*Vipv*(Ia + Is)/para.N;
        dVopv = -para.k*para.beta*Vopv*(Ia + Is)/para.N - para.delta*Vopv;
        dStil = para.gamma*(Ia + Is) - para.k*para.beta*Stil*(Ia + Is)/para.N;
        %Equation for IsCount is same as for Is but we remove recoveries
        dIsCount = para.a*para.beta*S*(Ia + Is)/para.N;

        dPop = [dS; dIs; dIa; dVipv; dVopv; dStil; dIsCount];

    end

end
