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

nstep  =  4000;
nprint =    50;
dtime  = 0.2;        

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
%
pix=4.0*atan(1.0);

%%  Create initial level set function for initial void
cc=[Nx/2 Ny/2];
R=50;
% cc=[3 3; 8 3; 13 3;3 10;9 10; 4.5 17];
% v_rad=[1.7 1.7; 1.7 1.7;1.7 1.7;1.7 1.7;1.7 1.7;1.7 1.7];
% cc=[10 4; 2.5 11]; 
% v_rad=[2 1.3; 1.3 2.3];

% cc=[3 3; 8 3; 13 3;3 10;9 10; 4.5 17];
% v_rad=[1.7 1.7; 1.7 1.7;1.7 1.7;1.7 1.7;1.7 1.7;1.7 1.7];
%v_rad=[1 1,1 1; 1 1];
[X Y]= meshgrid(1:Nx,1:Ny);

phi =((X-cc(1)).^2 +(Y-cc(2)).^2)-R^2;     


contourf(phi,[0 0])
Vn0=-ones(Nx,Ny);

for istep =1:nstep

phiold =phi;

%---
% calculate the laplacians and epsilon:
%---

%--gradients of phi:

[phidx,phidy]=gradient(phi);

%Vn0=-0.01*phidx;
s=sqrt(phidx.^2 + phidy.^2);
phi=phiold + dtime*Vn0.*s;

if(mod(istep,10)==0)
    fprintf('done step: %5d\n',istep);
    phi=phi + dtime*sign(phi).*(1-s);
end
%--- write vtk file:

%write_vtk_grid_values(Nx,Ny,dx,dy,istep,phi,tempr);  
figure(1)
contourf(phi,[0 0])
end %istep

%--- calculate compute time:

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);
   
