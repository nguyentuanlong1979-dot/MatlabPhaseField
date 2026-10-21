%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase-field for C3S dissolution and CSH, CH precipitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%== get initial wall time:
clear all; close all;
warning('off', 'MATLAB:contour:ConstantData');
colorblue= [0.06,0.06,1];
time0=clock();
format long;
%-- Simulation cell parameters:
Nx = 300;
Ny = 400;
NxNy= Nx*Ny;
dx = 0.03;
dy = 0.03; 
%--- Time integration parameters:
nstep  =  80000;
nprint =    200;
nnuc =    20;
dtime= 2.e-4;        
Time=0:dtime:(nstep-1)*dtime;
pix=4.0*atan(1.0);
% Gamma=0.05;

%--- Material specific parameters:
% phase 1: C3S
% phase 2: CSH
% phase 3: water - only for Temperature and Concentration diffusion
Lsurf = [1.0e+2 1.0e+2];
Alpha_C = [0.002 .015 .18];
% Alpha_T=[0.4 0.6 1.2];
Alpha_T=[0.01 0.02 0.04];

% In a PhD thesis summary (Sowoidnich) the authors list reaction enthalpies for alite 
% hydration steps: C₃S dissolution highly exothermic (≈ −132 kJ/mol), 
% C–S-H formation reported as endothermic, and portlandite (CH) precipitation exothermic (+25 kJ/mol)
Kappa_T = [-2.8 1.5];
% Kappa_T = [-3.8 1.2];

Kappa_C=[-0.81 -0.55];
Kappa_Int=[0.1 5.0];

epsilonb= [0.02 0.015];
theta0 = [0.00 0.00];  
delta = [0.05 0.005]; % strength of anisotropic
aniso = [0.1 0.1]; 
Ea=[35000 42000]; % J/mol %an icrease of temperature leads to increase of dissolution and decrase of precipitation
A=[1.0e6 2e7]; % called the frequency factor:
Rg=8.314; %J/mol K  % This is relate to Activation energy Ea  
kv=273.15; %Kevin number
a=0.1;
%Teq = [21.0 21.0];
Ceq=[24.0 20.0];
ConSolid=24.01;
tempr = 19.5*ones(Nx,Ny);
con = 19.9*ones(Nx,Ny);

%--- Initiliaze and introduce initial nuclai:
C3SThick=50;
SeedPoint=15.0;
Sradius=2.0;
Nphase=2;
RateSeed=400;
StartSeed=2500;

[phi, PointL, PointU]=Initial_StatePlainP(Nphase,Nx,Ny,C3SThick,SeedPoint);
inrange =(phi(:,:,1) >= 0.99);
con(inrange)=ConSolid;

index=0;
% LboundL=C3SThick-Sradius;
% UboundL=C3SThick+Sradius;

% phi(C3SThick:Rbound,Point(index)-Sradius:Point(index)+Sradius,2)=0.5;  
% CSH=phi(:,:,2);

% [phi]=Seeding(Nphase,Nx,Ny,C3SThick,Sradius,SeedPoint);

%save('InitialPhi.mat','C3S','phi');
% Daten = matfile('InitialPhi.mat');

%--- Laplacian templet:
% [laplacian1] =laplacianA(Nx,Ny,dx,dy);
% [~,~,laplacian1] = laplacian([Nx Ny], {'NN', 'P'});
[~,~,laplacian1] = laplacian([Nx Ny], {'P', 'P'});

% [~,~,laplacian1] = laplacian([Nx Ny]);

laplacian1=-laplacian1/(dx*dy);
%--- Evolution
C3Sini=phi(:,:,1);

Vol_phase= zeros(nstep,Nphase);
Rate_phase= zeros(nstep,Nphase);
% SeedPoint=1.0;
% TimeRate = zeros(nstep/5,1);
% Rate_phaseN = zeros(nstep/5,Nphase);
Calorie=zeros(nstep,1);

for istep =1:nstep
    %inrange =(C3S >= 0.9);
%    if(mod(istep,nnuc) == 0 && istep<300 )
%       [Nuclear]=Initial_State2(Nphase,Nx,Ny,R,Sradius,SeedPoint);
%        phi(:,:,2)= phi(:,:,2) + Nuclear(:,:,2);
%        phi(phi(:,:,2)>1)=1.0;
%    end

    LatentHeat=0;
    solid=0;
    Transform=0; %additive term
    for n=1:Nphase
        phin=phi(:,:,n);
        phinOld=phin;
        % calculate the laplacians and epsilon:
        lap_phi2 = laplacian1*reshape(phin,NxNy,1);
        [lap_phi]=reshape(lap_phi2,Nx,Ny);
        %--gradients of phi:
        [phidy,phidx]=gradient_matD(phin,Nx,Ny,dx,dy);
        %-- calculate angle:
        % theta =atan2(phidy,phidx);
        % theta =0;
%         Noise=(1+a*(rand(Nx,Ny)-0.5));
        %epsilonN = epsilonb(n).*(1+a*(rand(Nx,Ny)-0.5));
        % epsilon = epsilonb(n).*(1.0+delta(n)*cos(aniso(n)*(theta-theta0(n))));
        % epsilon_deriv = -aniso(n)*delta(n)*epsilonb(n).*sin(aniso(n).*(theta-theta0(n)));
        
        %--- first term:
        % dummyx =epsilon.*epsilon_deriv.*phidx;
        % [term1,dummy] =gradient_matD(dummyx,Nx,Ny,dx,dy);
        
        %--- second term:
        % dummyy =-epsilon.*epsilon_deriv.*phidy;
        % [dummy,term2] =gradient_matD(dummyy,Nx,Ny,dx,dy);
        
        %--- factor m_alpha:
        k_C= A(n).*exp(-Ea(n)./Rg./(tempr+kv));   
%         end
        m_alpha =(1/pix)*atan(k_C(n).*(con-Ceq(n)));
        % m_alpha =k_C(n).*(con-Ceq(n));
        % m_alpha(m_alpha>0.49) = 0.499;
        % m_alpha(m_alpha<-0.49) = -0.499;

        [dfdphi] = NonEq_Free_Energy(Nx,Ny,phin,phi,m_alpha,n,Nphase,Kappa_Int(n),a);
        
        % DriveForce=term1 +term2 + epsilon.^2 .* lap_phi - dfdphi;
        DriveForce= epsilonb(n).^2 .* lap_phi - dfdphi;
        

        %DriveForce(inrange)=0.0;
        phin =phin +(dtime*Lsurf(n))*DriveForce;       
        solid=solid + phin;
        LatentHeat=LatentHeat + Kappa_T(n)*(phin-phinOld);
        Transform=Transform + Kappa_C(n)*(phin-phinOld);
        Vol_phase(istep,n)= sum(sum(phin));
        Rate_phase(istep,n)=sum(sum(phin-phinOld));
    % %boundary Conditions No flux boundary 
    % phin=Boundary_Neumann(phin,Nx,Ny);
    phi(:,:,n)=phin;
    end
    
    %-- evolve temperature in medium
    temp2 = reshape(tempr,NxNy,1);
    lap_temp2=laplacian1*temp2;
    [lap_tempr] =reshape(lap_temp2,Nx,Ny);
    
    water=(solid<0.99);
    alpha_total=Alpha_T(1)*phi(:,:,1) + Alpha_T(2)*phi(:,:,2) + Alpha_T(3)*water;
    MeanTempOld= mean(tempr,'all');
    tempr = tempr + dtime*alpha_total.*lap_tempr + LatentHeat;
    Calorie(istep)= mean(tempr,'all') -  MeanTempOld;
    
    %---evolution of concentration
    con1 = reshape(con,NxNy,1);
    
    lap_con1=laplacian1*con1;
    [lap_con] =reshape(lap_con1,Nx,Ny);
        % 
        % if(mod(istep-1,5) == 0 ) && istep~=1
        % TimeRate((istep-1)/5)=Time(istep);
        % Rate_phaseN((istep-1)/5,n)=((Rate_phase(istep,n) - Rate_phase(istep-5,n)))/(5*dtime);
        % end

    alpha_total=Alpha_C(1)*phi(:,:,1) + Alpha_C(2)*phi(:,:,2) + Alpha_C(3)*water;
    con = con + dtime*alpha_total.*lap_con + Transform;
   
    % %boundary Conditions No flux boundary 
    % tempr=Boundary_Neumann(tempr,Nx,Ny);
    % con=Boundary_Neumann(con,Nx,Ny);
    if istep>StartSeed
    if(mod(istep,RateSeed) == 0 ) && index < length(PointL)
        index=index+1;
        phi(C3SThick-Sradius:C3SThick+Sradius,PointL(index)-Sradius:PointL(index)+Sradius,2)=0.8;  
        phi(Nx-C3SThick-Sradius:Nx-C3SThick+Sradius,PointU(index)-Sradius:PointU(index)+Sradius,2)=0.8;  
    end
    end
%---- print results
    if(mod(istep-1,nprint) == 0 )
        fprintf('done step: %5d\n',istep);
        %--- write vtk file:
        %wite_vtk_grid_values(Nx,Ny,dx,dy,istep,phi,tempr);
        figure(1)
        contourf(phi(:,:,1),1, 'LineColor', colorblue);
        hold on
        contourf(phi(:,:,2), [.99 .99], 'r');
        hold on
        %contour(C3Sini, 'LineColor', 'm','LineWidth',0.2); 
        yline([C3SThick Nx-C3SThick], '--r', 'LineWidth',.2); 
        axis equal
        %contour(GridX,GridY,-reshape(PhiIni,ex+1,ey+1)',[0 0],,'LineWidth',1);
        %contour(GridX,GridY,-reshape(PhiIni,ex+1,ey+1)',[0 0],,'LineWidth',1);
        FigTitle=strcat('Phase time :  ',num2str(istep*dtime));
        title(FigTitle)
        fignum=strcat('Figs/Phase/',num2str((istep-1)/nprint),'.png');      
        % F=getframe(gcf);
        % imwrite(F.cdata,fignum)
        % saveas(figure(1),fignum,'png')  
       
        figure(2)
        contourf(tempr);
        title('Temperature distribution')
        colorbar('southoutside')
        caxis([19.0 21])
        axis equal       
        FigTitle=strcat('Temp. time :  ',num2str(istep*dtime));
        title(FigTitle)
        fignum=strcat('Figs/Temp/',num2str((istep-1)/nprint),'.png');
        % F=getframe(gcf);
        % imwrite(F.cdata,fignum)
        % saveas(figure(2),fignum,'png')  
        
        figure(3)
        contourf(con);
        title('Concentration distribution')
        colorbar('southoutside')
        axis equal
        caxis([20 24])       
        FigTitle=strcat('Concent time :  ',num2str(istep*dtime));
        title(FigTitle)
        fignum=strcat('Figs/Concentration/',num2str((istep-1)/nprint),'.png');
        F=getframe(gcf);
        imwrite(F.cdata,fignum)
        % saveas(figure(3),fignum,'png')  
        % pause(0.01)
    end %if
    
end %istep

    figure(4)
    mesh(tempr); % for a better view, the LSF is displayed upside down
    hold on;  
    contour(tempr, 1, 'r','LineWidth',2);
    view([80 35]);
    axis on;
    [nrow, ncol]=size(tempr);
    axis([1 ncol 1 nrow 19 21]);
    set(gca,'ZTick',[19:1:21]);
    set(gca,'FontSize',14)

