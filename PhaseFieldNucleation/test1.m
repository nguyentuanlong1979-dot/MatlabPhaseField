clear all
Ea=[35000 42000]; % J/mol %an icrease of temperature leads to increase of dissolution and decrase of precipitation
A=[6.0e6 3e7]; % called the frequency factor:
Rg=8.314; %J/mol K  % This is relate to Activation energy Ea  
kv=273.15; %Kevin number
pix=4.0*atan(1.0);
Ceq=[24.0 20.0];

con= 15:0.1:25;
tempr=15:0.1:25;
% tempr=20;
% con= 21;

n=1;
 con_hat= (con - Ceq(n))./(Ceq(1) - Ceq(2));
k_C1= A(n).*exp(-Ea(n)./Rg./(tempr+kv));
m_alpha1 =(1/pix)*atan(k_C1.*con_hat);

n=2;
con_hat= (con - Ceq(n))./(Ceq(1) - Ceq(2));
k_C2= A(n).*exp(-Ea(n)./Rg./(tempr+kv));
m_alpha2 =(1/pix)*atan(k_C2.*con_hat);

figure
plot(tempr, k_C1, '-k')
hold on
plot(tempr, k_C2, '-r')


figure
plot(con, m_alpha1, '-k')
hold on
plot(con, m_alpha2, '-r')


xlim([20 24])

% figure
% t=-10:0.001:-1;
% E=exp(t);
% plot(t,E)