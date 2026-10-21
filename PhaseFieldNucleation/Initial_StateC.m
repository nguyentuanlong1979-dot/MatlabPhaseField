function [C3S,CSH]=Initial_StateC(Nx,Ny,R)
%== get intial wall time:
%-- Simulation cell parameters:
% Nx = 300;
% Ny = 300;
% R = 30.0;
% numPoints = 8;
% RCH=5;
% PosCH = [4*Nx/6 4*Ny/6];
xc1 = Nx/4;
yc1 = Ny/4;

xc2 = 3*Nx/4;
yc2 = 3*Ny/4;

%CSH  =zeros(Nx,Ny);
%C3S  =CSH;
%CH   =C3S;

[X,Y] = meshgrid(1:Nx,1:Ny);

xx1 =sqrt((X-xc1).^2 +(Y-yc1).^2);
C3S= (xx1<=R);
C3S(C3S>=1)=0.999;
%%
xx1 =sqrt((X-xc2).^2 +(Y-yc2).^2);

CSH= (xx1<=R);

CSH(CSH>=1)=0.999;

end
    