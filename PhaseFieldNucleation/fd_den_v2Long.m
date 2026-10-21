%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PHASE-FIELD FINITE-DIFFRENCE CODE FOR      %
%             DENDRITIC SOLIDIFICATION          %
%           (OPTIMIZED FOR MATLAB/OCTAVE)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab file tries to use the ordinary laplacian function
clear all; close all;
%== get initial wall time:
time0=clock();
format long;

%-- Simulation cell parameters:

Nx = 300;
Ny = 300;
NxNy= Nx*Ny;

dx = 0.03;
dy = 0.03; 
X=Nx*dx;
Y=Ny*dy;

%--- Time integration parameters:

nstep  =  8000;
nprint =    50;
dtime  = 1.e-6;   

%--- Material specific parameters:

tau = 0.0003; 
epsilonb = 0.01; 
mu = 1.0;
kappa = 1.8;    
delta = 0.005; 
aniso = 8.0; 
%aniso = randi([4 12],1,nstep);
alpha = 0.9; 
gamma = 10.0;
teq = 1.0;
theta0 = 0.2;  
seed = 5.0;
%
pix=4.0*atan(1.0);

%--- Initiliaze and introduce initial nuclai:

%[phi,tempr] = nucleus(Nx,Ny,seed);
phi = zeros(Nx,Ny);
tempr = zeros(Nx,Ny);

CenterPoint=Nx/2;
for i = 1:6
    for j = CenterPoint-3:CenterPoint+3
    phi(i,j) = 1; % Makes a portion of the domain on one side solid
    end   
end

% for j = 1:6
%     for i = CenterPoint-3:CenterPoint+3
%     phi(i,j) = 1; % Makes a portion of the domain on one side solid
%     end   
% end

%--- Laplacian templet:

%[laplacian1] =laplacian1(Nx,Ny,dx,dy);
[~,~,lapl] = laplacian([Nx Ny], {'DD', 'DD'});


%---
%--- Evolution
%---
% phi2 = reshape(phi,1,NxNy);
% lap_phi2 = laplacian(phi2,1);

for istep =1:nstep
phiold =phi;

%---
% calculate the laplacians and epsilon:
%---

phi2 = reshape(phi',NxNy,1);

lap_phi2= (lapl/(dx*dy))*phi2;

[lap_phi]=vec2matx(lap_phi2,Nx);

%--

tempx = reshape(tempr',NxNy,1);
% 
lap_tempx =(lapl/(dx*dy))*tempx;
% 
[lap_tempr] =vec2matx(lap_tempx,Nx);

%--gradients of phi:

[phidx,phidy]=gradient(phi,dx,dy);

%[phidy1,phidx1]=gradient_mat(phi,Nx,Ny,dx,dy);


%-- calculate angle:

theta =atan2(phidy,phidx);

%--- epsilon and its derivative:

epsilon = epsilonb*(1.0+delta*cos(aniso*(theta-theta0)));

epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(theta-theta0));

%--- first term:

dummyx =epsilon.*epsilon_deriv.*phidx;

[dummy,term1] =gradient(dummyx,dx,dy);

%--- second term:

dummyy =-epsilon.*epsilon_deriv.*phidy;

[term2,dummy] =gradient(dummyy,dx,dy);

%--- factor m:

m =(alpha/pix)*atan(gamma*(teq-tempr));

%-- Time integration:

phi = phi +(dtime/tau) *(term1 +term2 + epsilon.^2 .* lap_phi + ...
			 phiold.*(1.0-phiold).*(phiold - 0.5 + m));

phi(1:Nx,1)=0;
phi(1,1:Ny)=0;
phi(1:Nx,Ny)=0;
phi(Nx,1:Ny)=0;
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
contourf(phi)

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

    figure(3);
    mesh(phi); % for a better view, the LSF is displayed upside down
    hold on;  contour(phi, [0,0], 'r','LineWidth',2);
    view([-80 35]);
    axis on;
    [nrow, ncol]=size(phi);
    axis([1 ncol 1 nrow -5 5]);
    set(gca,'ZTick',[-3:1:3]);
    set(gca,'FontSize',14)
    
