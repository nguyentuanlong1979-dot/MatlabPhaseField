function [phi]=Seeding(Nphase,Nx,Ny,C3SThick,Dcsh,numPoints)
%== get intial wall time:
%-- Simulation cell parameters:
% clear all; close all;
% Nphase=2;
% Nx = 300;
% Ny = 300;
% C3SThick = 30;
% numPoints = 8;
% Dcsh=4;

Rbound=C3SThick+2*Dcsh;

phi = zeros(Nx,Ny,Nphase);

% C3S 
phi(1:C3SThick,:,1)=0.999;


Center = ceil(linspace(0,Ny,2*numPoints)); % 360 is the total number of points
Center(1)=[];
Center(end)=[];

randomIndexes = randperm(length(Center), numPoints);

 for i=1:numPoints
 Point=Center(randomIndexes(i));
 phi(C3SThick:Rbound,Point-Dcsh:Point+Dcsh,2)=0.8;  
 end

% CSH=phi(:,:,2);
% C3S=phi(:,:,1);
% contourf(C3S)
% hold on
% contourf(CSH,[.99 .99])
end

