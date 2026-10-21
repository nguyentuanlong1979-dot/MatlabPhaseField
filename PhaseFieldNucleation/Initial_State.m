function [C3S,CSH,CH]=Initial_State(Nx,Ny,R,Rch,numPoints)
%== get intial wall time:
%-- Simulation cell parameters:
% Nx = 300;
% Ny = 300;
% R = 30.0;
% numPoints = 8;

xc = Nx/2;
yc = Ny/2;
Rbound=R+3;
Sradius=3;
xch = 4*Nx/6;
ych = 4*Ny/6;

CSH  =zeros(Nx,Ny);
C3S  =CSH;
CH   =C3S;

for i=1:Nx
    for j=1:Ny
        xx1 =sqrt((i-xc)^2 +(j-yc)^2);
        xx2 =sqrt((i-xch)^2 +(j-ych)^2);
        if( xx1 <= R)
            C3S(i,j) = 0.999;
        end
        if( xx2 <= Rch)
            CH(i,j) = 0.999;
        end
    end %if
end % j

[xcr, ycr] = GetRandomPointsOnACircle(xc, yc, Rbound, numPoints);
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
    