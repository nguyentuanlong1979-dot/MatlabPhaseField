function [dfdphi] = NonEq_Free_EnergyCEM(Nx,Ny,phii,phi,phase,Nphase,kappa,W,c,T,G0)
format long;
sum_phi =zeros(Nx,Ny);
Rg=8.314; %J/mol K  % This is relate to Activation energy Ea  
kv=273.15; %Kevin number
Ceq=[25.0 25.0];
c=c./Ceq(phase);
T=(T+kv);
for iphase=1:Nphase
    if(iphase ~= phase)
    sum_phi = sum_phi + phi( :, :,iphase).^2;
    end
end

% Noise=(1+a*(rand(Nx,Ny)-0.5));

dgphi= phii.^3 - 3/2*phii.^2 + 1/2*phii;
dpphi=30*(phii.^2).*(1-phii).^2;
Gi0 = (1-c).*G0 + (20500-9.86*T).*c.*(1-c) + Rg*T.*((c.^2).*log(c) + (1-c).*log(1-c));
if phase==1
  Gi0=-Gi0;
else
    Gi0=Gi0;
end
%Gi0=-Rg*T*
% dfdphi = phii.*(phii-1).*(phii - 0.5 + m) + 2*kappa*phii.*sum_phi;

dfdphi = W*T.*dgphi - Gi0.*dpphi + 2*kappa*phii.*sum_phi;

% dfdphi=dfdphi./max(max(abs(dfdphi)));

end