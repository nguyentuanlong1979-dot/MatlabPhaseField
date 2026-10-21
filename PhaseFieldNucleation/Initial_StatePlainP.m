function [phi,PointL,PointU]=Initial_StatePlainP(Nphase,Nx,Ny,C3SThick,numPoints)
%== get intial wall time:
%-- Simulation cell parameters:
% clear all; close all;
% Nphase=2;
% Nx = 300;
% Ny = 300;
% C3SThick = 30;
% numPoints = 8;
% Dcsh=4;

phi = zeros(Nx,Ny,Nphase);

% C3S 
phi([1:C3SThick,Nx-C3SThick:Nx],:,1)=1.0;

Center = ceil(linspace(0,Ny,3*numPoints)); % 360 is the total number of points
Center(1)=[];
Center(end)=[];

randomIndexesU = randperm(length(Center), numPoints);
PointU=zeros(numPoints,1);
randomIndexesL = randperm(length(Center), numPoints);
PointL=zeros(numPoints,1);

 for i=1:numPoints
 PointL(i,1)=Center(randomIndexesL(i));
 PointU(i,1)=Center(randomIndexesU(i));
 % phi(C3SThick:Rbound,Point-Dcsh:Point+Dcsh,2)=0.8;  
 end

% CSH=phi(:,:,2);
% C3S=phi(:,:,1);
% contourf(C3S)
% hold on
% contourf(CSH,[.99 .99])
end

