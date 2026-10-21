%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%     PHASE-FIELD FINITE-DIFFERENCE CODE FOR        %
%              DENDRITIC SOLIDIFICATION             % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%== get intial wall time:
time0=clock();
format long;

%-- Simulation cell parameters:

Nx = 300;
Ny = 300;
NxNy= Nx*Ny;

dx = 0.03;
dy = 0.03; 

%--- Time integration parameters:

nstep  = 4000;
nprint =   50;
dtime  = 1.0e-4;   

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
seed = 5.0;
%
pix=4.0*atan(1.0);

%--Initiliaze and introduce initial nuclai:

[phi,tempr] = nucleus(Nx,Ny,seed);

%---
%- Evolution
%---

for istep =1:nstep

%----
% calculate the laplacians and epsilon:
%---

for i=1:Nx
for j=1:Ny

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

if(im == 0)
  im=Nx;
end
if(ip == (Nx+1))
  ip=1;
end

if(jm == 0) 
  jm = Ny;
end

if(jp == (Ny+1))
  jp=1;
end

hne=phi(ip,j);
hnw=phi(im,j);
hns=phi(i,jm);
hnn=phi(i,jp);
hnc=phi(i,j);

lap_phi(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);

hne=tempr(ip,j);
hnw=tempr(im,j);
hns=tempr(i,jm);
hnn=tempr(i,jp);
hnc=tempr(i,j);

lap_tempr(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);

%--gradients:

phidx(i,j) = (phi(ip,j) - phi(im,j))/(2.0*dx);
phidy(i,j) = (phi(i,jp) - phi(i,jm))/(2.0*dy);

%-- calculate angle:

theta = atan2(phidy(i,j),phidx(i,j));

%--- epsilon and its derivative:

epsilon(i,j) = epsilonb*(1.0+delta*cos(aniso*(theta-theta0)));

epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin(aniso*(theta-theta0));

end %for j
end %for i 

%----

for i=1:Nx
for j=1:Ny

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

jp=j+1;
jm=j-1;

ip=i+1;
im=i-1;

if(im == 0)
  im=Nx;
end
if(ip == (Nx+1))
  ip=1;
end

if(jm == 0) 
  jm = Ny;
end

if(jp == (Ny+1))
  jp=1;
end

phiold = phi(i,j);

%-- first term:

term1=   (epsilon(i,jp)*epsilon_deriv(i,jp)*phidx(i,jp) - ...
	  epsilon(i,jm)*epsilon_deriv(i,jm)*phidx(i,jm))/(2.0*dy);

%-- second term:

term2=  -(epsilon(ip,j)*epsilon_deriv(ip,j)*phidy(ip,j) - ... 
	  epsilon(im,j)*epsilon_deriv(im,j)*phidy(im,j))/(2.0*dx);

%-- factor m:

m = alpha/pi * atan(gamma*(teq-tempr(i,j)));

%-- Time integration:

phi(i,j) =phi(i,j) +(dtime/tau)*(term1 + term2 + epsilon(i,j)^2*lap_phi(i,j) + ...
				 phiold*(1.0-phiold)*(phiold -0.5 + m));

%-- evolve temperature:

tempr(i,j) =tempr(i,j) +dtime*lap_tempr(i,j) + kappa*(phi(i,j)-phiold);

end
end

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

write_vtk_grid_values(Nx,Ny,dx,dy,istep,phi,tempr);  

figure(1)
contourf(phi,1)
pause(0.01)
end %if

end %istep

%--- calculate compute time:

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);
   
