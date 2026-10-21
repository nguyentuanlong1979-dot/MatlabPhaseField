function [dfdphi] = NonEq_Free_EnergyN(Nx,Ny,phii,phi,m,phase,Nphase,kappa)
format long;
sum_phi =zeros(Nx,Ny);
X1 = rand(Nx,Ny) - 0.5; 
a = 0.02; % Amplitude of noise term
for iphase=1:Nphase
    if(iphase ~= phase)
    sum_phi = sum_phi + phi( :, :,iphase).^2;
    end
end

dfdphi = phii.*(1.0-phii).*(phii - 0.5 + m) - 2*kappa*phii.*sum_phi;

add_noise=a*phii*(1-phii).*X1;
inrange=(phii>0) - (phii>0.5);
dfdphi=dfdphi+add_noise.*inrange;

end