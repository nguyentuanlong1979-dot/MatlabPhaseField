
clear all; close all;
%== get intial wall time:
%-- Simulation cell parameters:
Nx = 300;
Ny = 300;
R = 30.0;
xc = Nx/2;
yc = Ny/2;

C3S  =zeros(Nx,Ny);


for i=1:Nx
for j=1:Ny

xx1 =sqrt((i-xc)^2 +(j-yc)^2);

if( xx1 <= R)
C3S(i,j) = 0.999;  
end

end %if

end % j

figure(1)
contourf(C3S)
% Random point on circle

numPoints = 8;
Rbound=R+3;
Sradius=3;
[xcr, ycr] = GetRandomPointsOnACircle(xc, yc, Rbound, numPoints);
CSH  =zeros(Nx,Ny);
for k=1:length(xcr)
    for i=1:Nx
        for j=1:Ny          
            xx1 =sqrt((i-xcr(k))^2 +(j-ycr(k))^2);     
            if( xx1 <= Sradius)
                CSH(i,j) = 0.999;
            end            
        end %if       
    end % j
end

phi  =zeros(Nx,Ny);

phi= CSH + C3S;

figure(4)
contourf(phi)



function [xRandom, yRandom] = GetRandomPointsOnACircle(xCenter, yCenter, radius, numPoints)
angles = linspace(0, 2*pi, 720); % 720 is the total number of points
x = radius * cos(angles) + xCenter;
y = radius * sin(angles) + yCenter;
% Now get random locations along the circle.
randomIndexes = randperm(length(angles), numPoints);
xRandom = x(randomIndexes);
yRandom = y(randomIndexes);
end

    