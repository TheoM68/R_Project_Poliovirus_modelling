% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provide a plot of Beckton's cases of both Sabin-like and VDPV cases in
% sewage during 2022
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Vector containing Sabin-like cases
SL = [0,0;14,0;28,2;42,0;56,0;70,1;91,5;105,5;119,8;126,7;133,3;140,5;146,7;147,8;154,5];
cVDPV = [0,0;14,0;28,0;42,0;56,0;70,0;91,0;105,0;119,0;126,0;133,1;140,3;146,0;147,0;154,2];

figure(1)
plot(SL(:,1),SL(:,2),'linewidth',2)
hold on
plot(cVDPV(:,1),cVDPV(:,2),'linewidth',2)
xlabel('Days since first recording')
ylabel('Isolates detected')
%title('Vaccine-like and cVDPV isolates detected in sewage data in Beckton recorded from 11 Jan 2022')
legend('Vaccine-Like','cVDPV','location','north west')