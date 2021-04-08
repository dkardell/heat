% File: conversion_line1AB.m
% Author: Dominik A. Kardell
% Date: December 2020

% Description: Estimates physical properties from FWI seismic velocity model using
% empirical relationships found in the literature. Writes out mat file
% required by hydrothermal modeling script.

%% Read binary file and feed parameters
clear, close all

fileID = fopen('line1AB_final.dir');
Vp = fread(fileID, [100 701], 'float');
fclose(fileID);
dx = 0.0375; %x-spacing in km
dy = 0.0375; %z-spacing in km
datum = 4.25; %Water depth at top in km
cutoff_vel = 3.1; %Sediment-crust cutoff velocity

%% Prepare plot parameters
xtic = 0:40:size(Vp,2); %Tick mark spacing for plots
ytic = 0:20:size(Vp,1);
xkm = xtic*dx; %Convert grid points to km
ykm = fliplr(ytic*dy)+datum;

%% Make depth grids for depth-dependent properties
%Depth below seafloor
z = zeros(size(Vp)); %Initialize depth matrix
for i = 1:size(Vp,1)
    z(i,:) = datum+i*dy;
end

%Sediment depth
zs = zeros(size(Vp)); %Initialize sed depth matrix
for i = 1:size(Vp,2)
    k = 0;
    for j = 1:size(Vp,1)
        if Vp(j,i) > 1.51
            zs(j,i) = k*dy;
            k = k+1;
        end
    end
end

%% Calculate porosity
%Using effective medium theory for crust
por_eff = zeros(size(Vp)); %Initialize porosity matrix
por_eff(Vp<1.51) = 1; %Water
por_eff(Vp>=1.51&Vp<cutoff_vel) = 0.72 - 0.987*zs(Vp>=1.51&Vp<cutoff_vel) + 0.830*zs(Vp>=1.51&Vp<cutoff_vel).^2; %Calcareous sediment (Hamilton, 1976)
% por_eff(Vp>=1.51&Vp<cutoff_vel) = 0.9 - 0.016*zs(Vp>=1.51&Vp<cutoff_vel) - 3.854*zs(Vp>=1.51&Vp<cutoff_vel).^2; %Radiolarian ooze (Hamilton, 1976)
% por_eff(Vp>=1.51&Vp<cutoff_vel) = 0.814 - 0.813*zs(Vp>=1.51&Vp<cutoff_vel) - 0.164*zs(Vp>=1.51&Vp<cutoff_vel).^2; %Pelagic clay (Hamilton, 1976)
% por_eff(Vp>=1.51&Vp<cutoff_vel) = 0.861 - 0.549*zs(Vp>=1.51&Vp<cutoff_vel) + 0.492*zs(Vp>=1.51&Vp<cutoff_vel).^2; %Diatomaceous ooze (Hamilton, 1976)
% por_eff(Vp>=1.51&Vp<cutoff_vel) = 0.909*(z(Vp>=1.51&Vp<cutoff_vel).^-0.073); %Seds (Spinelli & Fisher, 2014)
por_eff(Vp>=cutoff_vel) = -0.004*(Vp(Vp>=cutoff_vel).^3) + 0.029*(Vp(Vp>=cutoff_vel).^2) - 0.16*Vp(Vp>=cutoff_vel) + 0.76; %Crust (Marjanovic et al., 2019) / ?

%Using empirical relationships (Carlson, 2014) for crust
a_l=-0.129; b_l=0.181; c_l=0.0258; %Coefficients for lavas
a_d=-0.129; b_d=0.181; c_d=0.0258; %Coefficients for dikes
por_emp = zeros(size(Vp)); %Initialize porosity matrix
por_emp(Vp<1.51) = 1; %Water
por_emp(Vp>=1.51&Vp<cutoff_vel) = 0.72 - 0.987*zs(Vp>=1.51&Vp<cutoff_vel) + 0.830*zs(Vp>=1.51&Vp<cutoff_vel).^2; %Calcareous sediment (Hamilton, 1976)
por_emp(Vp>=cutoff_vel&Vp<5.16) = -b_l+sqrt((b_l^2)-4*a_l*(c_l-(1./(Vp(Vp>=cutoff_vel&Vp<5.16).^2)))./(2*a_l)); %Lavas
por_emp(Vp>=5.16) = -b_d+sqrt((b_d^2)-4*a_d*(c_d-(1./(Vp(Vp>=5.16).^2)))./(2*a_d)); %Dikes

%Difference between these two methods
por_diff = por_eff - por_emp;

%% Calculate permeability
K = zeros(size(Vp)); %Initialize permeability matrix
K(Vp < 1.51) = 0; %Water
K(Vp >= 1.51 & Vp < cutoff_vel) = (1.05e-18)*exp(2.17*(por_emp(Vp>=1.51&Vp<cutoff_vel)./(1-por_emp(Vp>=1.51&Vp<cutoff_vel)))); %Seds (Spinelli et al., 2004 hemipelagic)
K(Vp >= cutoff_vel & Vp < 5.16) = 10.^-(7.4+1.3*Vp(Vp>=cutoff_vel&Vp<5.16)); %Layer 2A (Carlson, 2014)
K(Vp >= 5.16) = 10.^(7-4*Vp(Vp>=5.16)); %Layer 2B (Carlson, 2014)

%% Calculate density
rho = zeros(size(Vp)); %Initialize density matrix
rho(Vp < 1.51) = 1.03; %Constant value for water

%Rock matrix density
rho(Vp >= 1.51 & Vp < cutoff_vel) = 2.71; %[g/cc] value for limestone from Hamilton (1885)
rho(Vp >= cutoff_vel) = 3.33; %[g/cc] value for basement from Hamilton (1885)

%Bulk density based on Vp
rhob = zeros(size(Vp)); %Initialize density matrix
rhob(Vp < 1.51) = 1.03; %Constant value for water
rhob(Vp >= 1.51 & Vp < cutoff_vel) = 1.512 + 1.631*zs(Vp>=1.51&Vp<cutoff_vel) - 1.373*zs(Vp>=1.51&Vp<cutoff_vel).^2; %Calcareous seds (Hamilton, 1976)
% rhob(Vp >= 1.51 & Vp < cutoff_vel) = 1.3813*(Vp(Vp>=1.51&Vp<cutoff_vel).^0.5083); %Adjusted Gardner eqn for shallow seds in Kerala-Konkan basin (Ojha & Sain, 2014)
rhob(Vp >= cutoff_vel) = 1.840*(Vp(Vp>=cutoff_vel).^0.25); %Adjusted Gardner eqn for mafic rocks (Brocher & Christensen, 2001)
%rhob = 0.31*((Vp*1000).^0.25); %Standard Gardner equation

%% Calculate thermal conductivity Ks
%For rock matrix (Hutchison, 1985)
Ks = zeros(size(Vp));
Ks(Vp < 1.51) = 0.6; %[1/(Wm*K)]
Ks(Vp >= 1.51 & Vp < cutoff_vel) = 2.93; %Value for limestone
Ks(Vp >= cutoff_vel) = 3.1;

%% Calculate specific heat capacity CPs
%For rock matrix (Hutchison, 1985)
CPs = zeros(size(Vp));
CPs(Vp < 1.51) = 3850; %[J/(kg*K)]
CPs(Vp >= 1.51 & Vp < cutoff_vel) = 1004; %Value for limestone
CPs(Vp >= cutoff_vel) = 1160;

%% Plot physical properties
figure(1)
subplot(4,2,1)
contourf(flipud(Vp), 100, 'LineStyle', 'none')
set(gca, 'colormap', jet); colorbar
caxis([1.5 7])
title('P-wave Velocity [km/s]')

subplot(4,2,3)
contourf(flipud(por_eff), 100, 'LineStyle', 'none')
jetf = flipud(jet);
set(gca, 'colormap', jetf); colorbar
caxis([0 0.7])
title('Porosity from Effective Medium Theory [frac]')

subplot(4,2,5)
contourf(flipud(por_emp), 100, 'LineStyle', 'none')
set(gca, 'colormap', jetf); colorbar
caxis([0 0.3])
title('Porosity from Carlson (2014) [frac]')

subplot(4,2,7)
contourf(flipud(por_diff), 100, 'LineStyle', 'none')
set(gca, 'colormap', jetf); colorbar
caxis([0 0.3])
title('Porosity Difference eff-emp [frac]')

subplot(4,2,2)
contourf(flipud(-log10(K)), 100, 'LineStyle', 'none')
logvec = fliplr(logspace(-10,-20,6));
set(gca, 'colormap', jet) 
caxis([10 20])
colorbar('YTick', -log10(fliplr(logvec)), 'YTickLabel', fliplr(logvec))
title('Permeability [m^2]')

subplot(4,2,4)
contourf(flipud(rho), 100, 'LineStyle', 'none')
set(gca, 'colormap', jet); colorbar
caxis([1 3])
title('Matrix Density [g/cm^3]')

subplot(4,2,6)
contourf(flipud(Ks), 100, 'LineStyle', 'none')
set(gca, 'colormap', jet); colorbar
%caxis([1 3])
title('Thermal Conductivity [1/(Wm*K)]')

subplot(4,2,8)
contourf(flipud(CPs), 100, 'LineStyle', 'none')
set(gca, 'colormap', jetf); colorbar
%caxis([1 3])
title('Specific Heat Capacity [J/(kg*K)]')

set(findobj(gcf,'type','axes'),'XTick',xtic,...
    'XTickLabel',xkm,'YTick',ytic,'YTickLabel',ykm)

%% Prepare and refine grids for flow modeling
wd = zeros(1,size(Vp,2)); %Initialize water depth vector
ocean_thick_vector = zeros(size(Vp,2),1);
for i = 1:size(Vp,2)
    k = 1;
    for j = 1:size(Vp,1)
        if Vp(j,i) > 1.51
            wd(i) = k;
            ocean_thick_vector(i) = datum*1000 + k*dx*1000;
            break
        end
        k = k+1;
    end
end

[K_wd,Ks_wd,rho_wd,por_emp_wd,CPs_wd,Vp_wd] = deal(zeros(size(Vp)));

for i = 1:length(wd)
    K_wd(1:end-wd(i)+1,i) = K(wd(i):end,i);
    Ks_wd(1:end-wd(i)+1,i) = Ks(wd(i):end,i);
    rho_wd(1:end-wd(i)+1,i) = rho(wd(i):end,i);
    por_emp_wd(1:end-wd(i)+1,i) = por_emp(wd(i):end,i);
    CPs_wd(1:end-wd(i)+1,i) = CPs(wd(i):end,i);
    Vp_wd(1:end-wd(i)+1,i) = Vp(wd(i):end,i);
end

K_wd = K_wd(1:end-max(wd),:);
Ks_wd = Ks_wd(1:end-max(wd),:);
rho_wd = rho_wd(1:end-max(wd),:);
por_emp_wd = por_emp_wd(1:end-max(wd),:); por_emp_wd(por_emp_wd <= 0) = 0.00001;
CPs_wd = CPs_wd(1:end-max(wd),:);
Vp_wd = Vp_wd(1:end-max(wd),:);

%% Export

n_div        = 0; % How many times to cut grid size in half
K_fine       = interp2(K_wd,n_div);
Ks_fine      = interp2(Ks_wd,n_div);
rho_fine     = interp2(rho_wd,n_div);
por_emp_fine = interp2(por_emp_wd,n_div);
CPs_fine     = interp2(CPs_wd,n_div);
Vp_fine      = interp2(Vp_wd,n_div);
% ocean_thick_vector = interp(ocean_thick_vector,2^n_div);

K_fine   = flipud(K_fine);
Ks_fine  = flipud(Ks_fine);
rho_fine = flipud(rho_fine);
por_fine = flipud(por_emp_fine);
CPs_fine = flipud(CPs_fine);
Vp_fine  = flipud(Vp_fine);

% save('phys_props.mat','rho_fine','por_fine','K_fine','Ks_fine','CPs_fine','Vp_fine','ocean_thick_vector','n_div');