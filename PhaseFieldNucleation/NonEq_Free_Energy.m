function [dfdphi] = NonEq_Free_Energy(Nx,Ny,phii,phi,m,phase,Nphase,kappa,a)
format long;
sum_phi =zeros(Nx,Ny);

for iphase=1:Nphase
    if(iphase ~= phase)
        sum_phi = sum_phi + phi( :, :,iphase).^2;
    end
end
    dfdphi = phii.*(phii-1).*(phii - 0.5 + m) + 2*kappa*phii.*sum_phi;

% Noise=(1+a*(rand(Nx,Ny)-0.5));
% if phase==2
%     dfdphi = phii.*(phii-1).*(phii - 0.5 + m) + 2*kappa*phii.*sum_phi;
% else
%     dfdphi = phii.*(phii-1).*(phii - 0.5 + m);
% end
% 
end