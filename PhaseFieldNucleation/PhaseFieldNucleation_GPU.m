%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPU-enabled phase-field for C3S dissolution and CSH/CH precipitation
% Based on the original PhaseFieldNucleation.m
%
% Requirements:
%   - MATLAB with gpuArray support
%   - Parallel Computing Toolbox for GPU computing
%   - Supporting functions on the MATLAB path:
%       Initial_StatePlainP.m
%       NonEq_Free_Energy.m
%
% Important:
%   The original sparse Laplacian multiplication has been replaced by a
%   periodic 5-point stencil using circshift. This matches the usual
%   2-D periodic finite-difference Laplacian, with the original dx*dy
%   scaling retained. Validate against the CPU version before production use.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
warning('off', 'MATLAB:contour:ConstantData');
colorblue = [0.06, 0.06, 1];
time0 = tic;
format long;

%% GPU setup
useGPU = true;   % Set false to run the array calculations on CPU
if useGPU
    try
        g = gpuDevice;
        fprintf('GPU enabled: %s\n', g.Name);
    catch ME
        warning('GPU initialization failed (%s). Falling back to CPU.', ME.message);
        useGPU = false;
    end
end

%% Simulation cell parameters
Nx = 300;
Ny = 400;
dx = 0.03;
dy = 0.03;

%% Time integration parameters
nstep  = 80000;
nprint = 200;
nnuc   = 20; %#ok<NASGU>
dtime  = 2.e-4;
Time   = 0:dtime:(nstep-1)*dtime; %#ok<NASGU>
pix    = 4.0*atan(1.0);

%% Material-specific parameters
% Phase 1: C3S; phase 2: CSH; phase 3: water for diffusion
Lsurf = [1.0e2 1.0e2];
Alpha_C = [0.002 .015 .18];
Alpha_T = [0.01 0.02 0.04];
Kappa_T = [-2.8 1.5];
Kappa_C = [-0.81 -0.55];
Kappa_Int = [0.1 5.0];
epsilonb = [0.02 0.015];
theta0 = [0.00 0.00]; %#ok<NASGU>
delta = [0.05 0.005]; %#ok<NASGU>
aniso = [0.1 0.1]; %#ok<NASGU>
Ea = [35000 42000]; % J/mol
A = [1.0e6 2e7];   % Frequency factors
Rg = 8.314;        % J/(mol K)
kv = 273.15;       % Celsius-to-Kelvin offset
a = 0.1;
Ceq = [24.0 20.0];
ConSolid = 24.01;

%% Initial conditions
tempr = 19.5*ones(Nx,Ny);
con   = 19.9*ones(Nx,Ny);

C3SThick = 50;
SeedPoint = 15.0;
Sradius = 2.0;
Nphase = 2;
RateSeed = 400;
StartSeed = 2500;

[phi, PointL, PointU] = Initial_StatePlainP( ...
    Nphase, Nx, Ny, C3SThick, SeedPoint);

inrange = (phi(:,:,1) >= 0.99);
con(inrange) = ConSolid;

index = 0;

%% Move large arrays to GPU
if useGPU
    phi   = gpuArray(phi);
    tempr = gpuArray(tempr);
    con   = gpuArray(con);
end

% Histories remain on GPU during the time loop to avoid repeated CPU/GPU
% synchronization from scalar history assignments.
if useGPU
    Vol_phase = gpuArray.zeros(nstep,Nphase);
    Rate_phase = gpuArray.zeros(nstep,Nphase);
    Calorie = gpuArray.zeros(nstep,1);
else
    Vol_phase = zeros(nstep,Nphase);
    Rate_phase = zeros(nstep,Nphase);
    Calorie = zeros(nstep,1);
end

%% Evolution
for istep = 1:nstep
    LatentHeat = 0;
    solid = 0;
    Transform = 0;

    for n = 1:Nphase
        phin = phi(:,:,n);
        phinOld = phin;

        % Periodic 5-point Laplacian (GPU-compatible).
        lap_phi = laplacianPeriodicGPU(phin, dx, dy);

        % The original gradient calculation is not used in the active
        % anisotropy terms, so it is intentionally omitted.

        k_C = A(n).*exp(-Ea(n)./Rg./(tempr+kv));
        m_alpha = (1/pix)*atan(k_C.*(con-Ceq(n)));

        % GPU-safe equivalent of NonEq_Free_Energy.m. Starting with a
        % scalar avoids mixing a CPU zeros(Nx,Ny) array with GPU arrays.
        sum_phi = 0;
        for iphase = 1:Nphase
            if iphase ~= n
                sum_phi = sum_phi + phi(:,:,iphase).^2;
            end
        end
        dfdphi = phin.*(phin-1).*(phin-0.5+m_alpha) ...
               + 2*Kappa_Int(n)*phin.*sum_phi;

        DriveForce = epsilonb(n).^2 .* lap_phi - dfdphi;
        phin = phin + (dtime*Lsurf(n))*DriveForce;

        solid = solid + phin;
        LatentHeat = LatentHeat + Kappa_T(n)*(phin-phinOld);
        Transform = Transform + Kappa_C(n)*(phin-phinOld);

        Vol_phase(istep,n) = sum(phin,'all');
        Rate_phase(istep,n) = sum(phin-phinOld,'all');
        phi(:,:,n) = phin;
    end

    %% Temperature evolution
    lap_tempr = laplacianPeriodicGPU(tempr, dx, dy);
    water = (solid < 0.99);
    alpha_total = Alpha_T(1)*phi(:,:,1) + ...
                  Alpha_T(2)*phi(:,:,2) + Alpha_T(3)*water;

    MeanTempOld = mean(tempr,'all');
    tempr = tempr + dtime*alpha_total.*lap_tempr + LatentHeat;
    Calorie(istep) = mean(tempr,'all') - MeanTempOld;

    %% Concentration evolution
    lap_con = laplacianPeriodicGPU(con, dx, dy);
    alpha_total = Alpha_C(1)*phi(:,:,1) + ...
                  Alpha_C(2)*phi(:,:,2) + Alpha_C(3)*water;
    con = con + dtime*alpha_total.*lap_con + Transform;

    %% Nucleation
    if istep > StartSeed
        if mod(istep,RateSeed) == 0 && index < length(PointL)
            index = index + 1;
            phi(C3SThick-Sradius:C3SThick+Sradius, ...
                PointL(index)-Sradius:PointL(index)+Sradius,2) = 0.8;
            phi(Nx-C3SThick-Sradius:Nx-C3SThick+Sradius, ...
                PointU(index)-Sradius:PointU(index)+Sradius,2) = 0.8;
        end
    end

    %% Plotting and image output: gather only at output intervals
    if mod(istep-1,nprint) == 0
        fprintf('done step: %5d\n', istep);

        phiPlot = gatherIfGPU(phi,useGPU);
        tempPlot = gatherIfGPU(tempr,useGPU);
        conPlot = gatherIfGPU(con,useGPU);

        figure(1);
        contourf(phiPlot(:,:,1),1,'LineColor',colorblue);
        hold on;
        contourf(phiPlot(:,:,2),[.99 .99],'r');
        yline([C3SThick Nx-C3SThick],'--r','LineWidth',.2);
        axis equal;
        title(strcat('Phase time : ',num2str(istep*dtime)));
        hold off;

        figure(2);
        contourf(tempPlot);
        title('Temperature distribution');
        colorbar('southoutside');
        caxis([19.0 21]);
        axis equal;
        title(strcat('Temp. time : ',num2str(istep*dtime)));

        figure(3);
        contourf(conPlot);
        title('Concentration distribution');
        colorbar('southoutside');
        axis equal;
        caxis([20 24]);
        title(strcat('Concent time : ',num2str(istep*dtime)));

        % Preserve the original concentration-image output behavior.
        outDir = fullfile('Figs','Concentration');
        if ~exist(outDir,'dir')
            mkdir(outDir);
        end
        F = getframe(gcf);
        imwrite(F.cdata, fullfile(outDir, ...
            [num2str((istep-1)/nprint) '.png']));
    end
end

%% Final plots
temprPlot = gatherIfGPU(tempr,useGPU);
figure(4);
mesh(temprPlot);
hold on;
contour(temprPlot,1,'r','LineWidth',2);
view([80 35]);
axis on;
[nrow,ncol] = size(temprPlot);
axis([1 ncol 1 nrow 19 21]);
set(gca,'ZTick',19:1:21);
set(gca,'FontSize',14);

% Bring histories back to CPU once, after the simulation.
Vol_phase = gatherIfGPU(Vol_phase,useGPU); %#ok<NASGU>
Rate_phase = gatherIfGPU(Rate_phase,useGPU); %#ok<NASGU>
Calorie = gatherIfGPU(Calorie,useGPU); %#ok<NASGU>

fprintf('Total elapsed time: %.3f seconds\n', toc(time0));

%% Local helper functions
function lap = laplacianPeriodicGPU(x, dx, dy)
%LAPLACIANPERIODICGPU Periodic 2-D five-point Laplacian.
% Uses circshift, which supports gpuArray inputs.
% Scaling follows the original code's dx*dy denominator.
    lap = (circshift(x,[1 0]) + circshift(x,[-1 0]) - 2*x)./(dx*dy) ...
        + (circshift(x,[0 1]) + circshift(x,[0 -1]) - 2*x)./(dx*dy);
end

function x = gatherIfGPU(x,useGPU)
    if useGPU
        x = gather(x);
    end
end
