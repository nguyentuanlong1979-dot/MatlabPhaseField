%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the effect of temperature in surface reaction rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all;
Ea=[-5000]; %J/mol
Rg=9.314; %J/mol K
A=[1.6];
A1=1/500;
T=20:.1:30; % celcius temperature
kv=273.15; %Kevin number

k=zeros(length(T),length(Ea));
for i=1:length(Ea)
k(:,i)= A(i).*exp(-Ea(i)./Rg./(T+kv));
end
k2= A1./Rg./(T+kv);

figure
plot(T,k)
xlabel('Temperature')
ylabel('k');

% figure
% plot(T,k2)
% xlabel('Temperature')
% ylabel('k');
