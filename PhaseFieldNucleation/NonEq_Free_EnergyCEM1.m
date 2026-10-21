function [dfdphi] = NonEq_Free_EnergyCEM1(Nx,Ny,phii,phi,phase,Nphase,kappa,W,c,T,G0)
format long;
sum_phi =zeros(Nx,Ny);
Rg=8.314; %J/mol K  % This is relate to Activation energy Ea  
kv=273.15; %Kevin number
xi=2.1e-2;
% Ceq=[24.0 20.0];
% c=(c-Ceq(1))./(Ceq(2)-Ceq(1));
T=(T+kv);
for iphase=1:Nphase
    if(iphase ~= phase)
    sum_phi = sum_phi + phi( :, :,iphase).^2;
    end
end

% Noise=(1+a*(rand(Nx,Ny)-0.5));

dgphi= phii.^3 - 3/2*phii.^2 + 1/2*phii;
dpphi=30*(phii.^2).*(1-phii).^2;
% Gi0 = (1-c).*G0 + (20500-9.86*T).*c.*(1-c) + Rg*T.*((c.^2).*log(c) + (1-c).*log(1-c));
% Gi0 = (1-c).*G0 + Rg*T.*(c.*log(c) + (1-c).*log(1 - c + 0.0111));

if phase==1
  Gi0=(1-c).*G0 + Rg*T.*(1-c).*log(1 - c);
else
  Gi0 = c.*G0 + Rg*T.*(c.*log(c));
end
%Gi0=-Rg
% *T*
% dfdphi = phii.*(phii-1).*(phii - 0.5 + m) + 2*kappa*phii.*sum_phi;

dfdphi = (W*T.*dgphi + Gi0.*dpphi + 2*kappa*phii.*sum_phi);
% dfdphi(abs(dfdphi)<1e-7)=0;
% AA=W*T.*dgphi;
% A= Gi0.*dpphi;

% dfdphi=dfdphi./max(max(abs(dfdphi)));

end