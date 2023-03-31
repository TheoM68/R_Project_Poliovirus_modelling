%Tau-leap model for our poliovirus model with quarantine of symptomatic
%infections
function [Classes] = Tauleap_polio_model_quarantine_Is(para,ICs,maxtime,tau)
    
    %Store starting point of simulation and
    %define starting time to be 0
    Classes = ICs;
    Classes.t = 0;
    %Introduce quarantine class.
    Classes.Q = 0; 
    %Introduce a counting class that counts symptomatic infections to
    %determine total cases of paralysis.
    Classes.IsCount = 0;
    %Introduce a counting class that counts how many times our population
    %becomes negative in our simulation and needs to be rectified.
    Classes.NegCount = 0;
    
    %Defin current state of model
    S = ICs.S;
    Is = ICs.Is;
    Ia = ICs.Ia;
    Vipv = ICs.Vipv;
    Vopv = ICs.Vopv;
    Stil = ICs.Stil;
    Q = 0;
    IsCount = 0;
    NegCount = 0;
    
    t = 0;
    
    %Run unitil maxtime is exceeded or there are no infections AND Vopv=0
    while ((t<maxtime) && ( (Ia>0)||(Is>0)||(Vopv>0) ))
        
        %Event rates
        r_1 = para.a*para.beta*S*(Ia + Is)/para.N; %New symtpomatic infection
        r_2 = (1-para.a)*para.beta*S*(Ia + Is)/para.N; %New asymptomatic infection (with no vaccine/infection history)
        r_3 = para.beta*Vipv*(Ia + Is)/para.N; %Infection after IPV vaccine
        r_4 = para.k*para.beta*Vopv*(Ia + Is)/para.N; %(Re-)infection after OPV vaccine
        r_5 = para.k*para.beta*Stil*(Ia + Is)/para.N; %Re-infection after infection with cVDPV2
        r_6 = para.nu*S; %Vaccination
        r_7 = para.delta*Vopv; %Mutation of live attenuated virus in OPV to cVDPV
        r_8 = para.gamma*Is; %Recovery from virus that induced paralytic poliomyelitis
        r_9 = para.gamma*Ia; %Recovery from virus which did not induce paralysis.
        %New events for quarantine
        r_10 = para.q*Is; %Quarantine of a symptomatic individual.
        r_11 = para.gamma*Q; %Recovery in quarantine - hence leave quarantine.
        
        %Compute how many events occur for each time step
        X_1 = poissrnd(r_1*tau);
        X_2 = poissrnd(r_2*tau);
        X_3 = poissrnd(r_3*tau);
        X_4 = poissrnd(r_4*tau);
        X_5 = poissrnd(r_5*tau);
        X_6 = poissrnd(r_6*tau);
        X_7 = poissrnd(r_7*tau);
        X_8 = poissrnd(r_8*tau);
        X_9 = poissrnd(r_9*tau);
        %New events for quarantine
        X_10 = poissrnd(r_10*tau);
        X_11 = poissrnd(r_11*tau);
        
        %Update classes
        S = S - X_1 - X_2 - X_6;
        Is = Is + X_1 - X_8 - X_10;
        Ia = Ia + X_2 + X_3 + X_4 + X_5 + X_7 - X_9;
        Vipv = Vipv - X_3 + X_6;
        Vopv = Vopv - X_4 - X_7;
        Q = Q + X_10 - X_11;
        Stil = Stil - X_5 + X_8 + X_9 + X_11;
        
        %Add all symptomatic infections to this new counting class
        IsCount = IsCount + X_1;
        
        %Check nothing is less than 0, if so "undo" one of the steps
        if S<0 %"undo" an infection OR vaccination
           NegCount = NegCount + 1;
           tmp = S;
           S = 0;
           
           %Randomly choose k of the negative susceptibles to be reverted
           %as infections, -tmp-k are reverted as vaccination
           k = binornd(-tmp, (X_1 + X_2)/(X_1 + X_2 + X_6));
           
           %Ensure k doesn't exceed the number of infections (X_1+X_2)
           k = min(k,X_1+X_2);
           %Ensure -tmp-k doesn't exceed the number of vaccinations (X_6)
           k = max(k,-tmp-X_6);
                      
           %Randomly choose m of the k negative susceptibles being reverted 
           %as infections to be reverted as symptomatic infections.
           %k-m will be reverted as asymptomatic infections.
           m = binornd(k,X_1/(X_1+X_2));
           
           %Ensure m doesn't exceed the number of symptomatic infections
           m = min(m,X_1);
           
           %Ensure k-m doesn't exceed the number asymptomatic infections
           m = max(m,k-X_2);
           
           %Update classes and revert infections and vaccinations
           Vipv = Vipv+tmp+k;
           Is = Is - m;
           Ia = Ia - (k-m);
           
        end
        
        if Is<0 %"undo" recovery
            NegCount = NegCount + 1;
            tmp = Is;
            Is = 0;
            Stil = Stil + tmp;
        
        end
        
        if Ia<0 %"undo" recovery
            NegCount = NegCount + 1;
            tmp = Ia;
            Ia = 0;
            Stil = Stil + tmp;
            
        end
        
        if Vipv<0 %"undo" infection
            NegCount = NegCount + 1;
            tmp = Vipv;
            Vipv = 0;
            Ia = Ia + tmp;
            
        end
        
        if Vopv<0 %"undo" infection OR mutation - both of these being 
                  %reverted decrease Ia so no random choice required.
            NegCount = NegCount + 1;
            tmp = Vopv;
            Vopv = 0;
            
            Ia = Ia + tmp;    
        
        end
        
        if Stil<0 %"undo" infection
            NegCount = NegCount + 1;
            tmp = Stil;
            Stil = 0;
            Ia = Ia + tmp;
            
        end
        
        if Q<0 %"undo" recovery in quarantine
            NegCount = NegCount + 1;
            tmp = Q;
            Q = 0;
            
            Stil = Stil + tmp;
        end
        
    %Update time
    t = t + tau;
    
    %Save information in the Classes structure
    Classes.t = [Classes.t t];
    Classes.S = [Classes.S S];
    Classes.Is = [Classes.Is Is];
    Classes.Ia = [Classes.Ia Ia];
    Classes.Vipv = [Classes.Vipv Vipv];
    Classes.Vopv = [Classes.Vopv Vopv];
    Classes.Stil = [Classes.Stil Stil];
    Classes.Q = [Classes.Q Q];
    Classes.IsCount = [Classes.IsCount IsCount];
    Classes.NegCount = [Classes.NegCount NegCount];
    end
        
end
        
        
        
        
        
        
        
        
        