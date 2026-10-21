function [C3S,CSH,CH]=Initial_State1(Nx,Ny,R,RCH,PosCH,numPoints)
%== get intial wall time:
%-- Simulation cell parameters:
% Nx = 300;
% Ny = 300;
% R = 30.0;
% numPoints = 8;
% RCH=5;
% PosCH = [4*Nx/6 4*Ny/6];
xc = Nx/2;
yc = Ny/2;
Rbound=R+3;
Sradius=3;

CSH  =zeros(Nx,Ny);
%C3S  =CSH;
%CH   =C3S;

L = linspace(0,2*pi,7);
xv = RCH*cos(L+pi()/15)'+PosCH(1);
yv = RCH*sin(L+pi()/15)'+PosCH(2);

[X,Y] = meshgrid(1:Nx,1:Ny);

[CH,~] = inpolygon(X,Y,xv,yv);

xx1 =sqrt((X-xc).^2 +(Y-yc).^2);

C3S= (xx1<=R);

C3S(C3S>=1)=0.999;

[xcr, ycr] = GetRandomPointsOnACircle(xc, yc, Rbound, numPoints);

for k=1:length(xcr)
CSHk=[];
xx1 =sqrt((X-xcr(k)).^2 +(Y-ycr(k)).^2);     
CSHk = (xx1 <= Sradius);
CSH=CSH + CSHk;
end
CSH(CSH>=1)=0.999;

function [xRandom, yRandom] = GetRandomPointsOnACircle(xCenter, yCenter, radius, numPoints)
angles = linspace(0, 2*pi, 720); % 720 is the total number of points
x = radius * cos(angles) + xCenter;
y = radius * sin(angles) + yCenter;
% Now get random locations along the circle.
randomIndexes = randperm(length(angles), numPoints);
xRandom = x(randomIndexes);
yRandom = y(randomIndexes);
end


end
    