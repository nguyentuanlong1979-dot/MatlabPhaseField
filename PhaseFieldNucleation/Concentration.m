clear all; close all;
colorblue= [0.06,0.06,1];

time0=clock();
format long;

%-- Simulation cell parameters:
Nx = 300;
Ny = 300;
NxNy= Nx*Ny;
dx = 0.03;
dy = 0.03; 

nstep  =  10000;
nprint =    50;
dtime  = 1.e-8;        
tau = 0.0003; 


%--- Initiliaze and introduce initial nuclai:
R=40;

%--- Initiliaze and introduce initial nuclai:
[C3S,CSH]=Initial_StateC(Nx,Ny,R);
 con = zeros(Nx,Ny);
 Diffusion=ones(Nx,Ny);
    figure(1)
    axis equal
    contourf(CSH, [0.99 0.99], 'LineColor', colorblue);
    hold on
    contourf(C3S, [0.99 0.99], 'r');

    [lap1] =laplacian1(Nx,Ny,dx,dy);
    
    [~,~,lap] = laplacian([Nx Ny], {'P', 'P'});
    lap=-lap/(dx*dy);

    BoundCSH=(CSH>0.99);
    BoundC3S=(C3S>0.99);
    
 for istep =1:nstep
    con(BoundC3S)=1;
    con(BoundCSH)=-1; 
    con = reshape(con,NxNy,1);   
    lap_con=lap*con;
    
    con= con +(dtime/tau)*lap_con;
    [con]=vec2matx(con,Nx);

    if(mod(istep,nprint) == 0 )
    fprintf('done step: %5d\n',istep);

    figure(2)
    axis equal
    contourf(con);
    title('Concentration distribution')
    colorbar('southoutside')
    pause(0.001)
    end
    
  end
    
    figure(3);
    mesh(con); % for a better view, the LSF is displayed upside down
    hold on;  
    contour(con, 1, 'r','LineWidth',2);
    view([80 35]);
    axis on;
    [nrow, ncol]=size(con);
    axis([1 ncol 1 nrow -2 2]);
    set(gca,'ZTick',-2:1:2);
    set(gca,'FontSize',14)