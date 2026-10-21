%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PHASE-FIELD FINITE-DIFFRENCE CODE FOR      %
%             DENDRITIC SOLIDIFICATION          %
%           (OPTIMIZED FOR MATLAB/OCTAVE)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%== get initial wall time:
clear all; close all;
time0=clock();
format long;

%-- Simulation cell parameters:

Nx = 300;
Ny = 300;
NxNy= Nx*Ny;

dx = 0.03;
dy = 0.03; 

%--- Time integration parameters:

nstep  =  15000;
nprint =    50;
dtime  = 1.e-4;        

%--- Material specific parameters:

tau = 0.0003; 
epsilonb = 0.01; 
mu = 1.0;
kappa = 1.8;    
delta = 0.02; 
aniso = 6.0; 
alpha = 0.9; 
gamma = 10.0;
teq = 1.0;
theta0 = 0.2;  
seed = 50.0;
a=0.001; 

%
pix=4.0*atan(1.0);

%--- Initiliaze and introduce initial nuclai:

[phi,tempr] = nucleus(Nx,Ny,seed);

%--- Laplacian templet:

[laplacian1] =laplacian1(Nx,Ny,dx,dy);

%---
%--- Evolution
%---

for istep =1:nstep

phiold =phi;

%---
% calculate the laplacians and epsilon:
%---
phi2 = reshape(phi',NxNy,1);

lap_phi2 = laplacian1*phi2;

[lap_phi]=vec2matx(lap_phi2,Nx);

%--

tempx = reshape(tempr',NxNy,1);

lap_tempx =laplacian1*tempx;

[lap_tempr] =vec2matx(lap_tempx,Nx);

%--gradients of phi:

[phidy,phidx]=gradient_mat(phi,Nx,Ny,dx,dy);

%-- calculate angle:

theta =atan2(phidy,phidx);

%--- epsilon and its derivative:

epsilon = epsilonb*(1.0+delta*cos(aniso*(theta-theta0)));

epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(theta-theta0));

%--- first term:

dummyx =epsilon.*epsilon_deriv.*phidx;

[term1,dummy] =gradient_mat(dummyx,Nx,Ny,dx,dy);

%--- second term:

dummyy =-epsilon.*epsilon_deriv.*phidy;

[dummy,term2] =gradient_mat(dummyy,Nx,Ny,dx,dy);

%--- factor m:

m =(alpha/pix)*atan(gamma*(teq-tempr));

%-- Time integration:
X1 = rand(Nx,Ny) - 0.5; 
phi = phi +(dtime/tau) *(term1 +term2 + epsilon.^2 .* lap_phi + ...
		    phiold.*(1.0-phiold).*(1.0-phiold));

%-- evolve temperature:
tempr = tempr + dtime*lap_tempr + kappa*(phi-phiold);


%---- print results

if(mod(istep,nprint) == 0 )

fprintf('done step: %5d\n',istep);

%fname1 =sprintf('time_%d.out',istep);
%out1 = fopen(fname1,'w');

%for i=1:Nx
%for j=1:Ny
%ii=(i-1)*Nx+j;
%fprintf(out1,'%5d %5d %14.6e %14.6e\n',i,j,phi(i,j),tempr(i,j))
%end
%end

%fclose(out1)

%--- write vtk file:

%write_vtk_grid_values(Nx,Ny,dx,dy,istep,phi,tempr);  
figure(1)
contourf(phi,1)

figure(2)
contourf(tempr)
title('Temperature distribution')
colorbar('southoutside')
pause(0.01)
end %if 

end %istep

%--- calculate compute time:

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);
   
