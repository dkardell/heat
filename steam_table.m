function [Rho, Cp, Mu] = steam_table(T,P,Grid)

%%
% author: Dominik A. Kardell
% modified from: Evan J. Ramos 
% date: December 2020

%Description: This function loads steam table data regarding fluid density,
%heat capacity, and dynamic viscosity and given a set of temperatures and
%pressures, returns interpolated values from matrices.

%Inputs:
%         T    --> Grid.N x 1 vector containing temperatures [K]
%         P    --> Grid.N x 1 vector containing pressures    [Pa]
%         Grid --> structure containing the grid properties

%Outputs:
%         Rho  --> Grid.N x 1 vector containing fluid density [kg/m^3]
%         Cp   --> Grid.N x 1 vector containing heat capacity [J/mol/K]
%         Mu   --> Grid.N x 1 vector containing viscosity     [Pa*s]

%% Declare persistent temperatures and enthalpies

% NOTE:
% By definition, heat capacity is the derivative of enthalpy with respect
% to temperature. Since this function will be used iteratively as
% temperature evolves, the temperatures and enthalpies from the previous
% time step will be stored in a persistent variable. Thus, the input
% pressures and temperatures (at that current time step) can be used to
% calculate the heat capacity at each point in space.
% 
persistent T_old
persistent H_old

if isempty(T_old)
    T_old = T;
    H_old = P;
end

%% Load steam tables

% NOTE:
% Density: from Pitzer et al. (1984)
% Enthalpy: from Burnham et al. (1969)
% Viscosity: from Haar et al. (1984)

load('density.mat');
enthalpy  = load('steam_table_enthalpy.txt');
viscosity = load('steam_table_viscosity.txt')';

%% Make interpolated matrices from steam table data

%Density
t_rho                    = density(2:end,1);
p_rho                    = density(1,2:end);
[P_rho, T_rho]           = meshgrid(p_rho,t_rho);

t_rho_fine               = linspace(min(t_rho),max(t_rho));
p_rho_fine               = linspace(min(p_rho),max(p_rho));
[P_rho_fine, T_rho_fine] = meshgrid(p_rho_fine,t_rho_fine);

rho_fine                 = interp2(P_rho,T_rho,density(2:end,2:end),...
                                   P_rho_fine,T_rho_fine,'spline');
rho_fine(rho_fine <= 400) = 400;

%Enthalpy
t_ent                    = enthalpy(2:end,1);
p_ent                    = enthalpy(1,2:end);
[P_ent, T_ent]           = meshgrid(p_ent,t_ent);

t_ent_fine               = linspace(min(t_ent),max(t_ent));
p_ent_fine               = linspace(min(p_ent),max(p_ent));
[P_ent_fine, T_ent_fine] = meshgrid(p_ent_fine,t_ent_fine);

ent_fine                 = interp2(P_ent,T_ent,enthalpy(2:end,2:end),...
                                   P_ent_fine,T_ent_fine,'spline');

%Viscosity
t_vis                    = viscosity(2:end,1);
p_vis                    = viscosity(1,2:end);
[P_vis, T_vis]           = meshgrid(p_vis,t_vis);

t_vis_fine               = linspace(min(t_vis),max(t_vis));
p_vis_fine               = linspace(min(p_vis),max(p_vis));
[P_vis_fine, T_vis_fine] = meshgrid(p_vis_fine,t_vis_fine);

vis_fine                 = interp2(P_vis,T_vis,viscosity(2:end,2:end),...
                                   P_vis_fine,T_vis_fine,'spline');

%% Locate temperature and pressure in fine vectors --> new fluid properties

Rho = zeros(Grid.N,1);
Cp  = zeros(Grid.N,1);
Mu  = zeros(Grid.N,1);

for i = 1:Grid.N
    %Density
    [~, r_rho] = min(abs(T(i)-t_rho_fine));
    [~, c_rho] = min(abs(P(i)-p_rho_fine));
    Rho(i)     = rho_fine(r_rho,c_rho);
    
% %     %Enthalpy
%     [~, r_ent] = min(abs(T(i)-t_ent_fine));
%     [~, c_ent] = min(abs(P(i)-p_ent_fine));
%     Ent        = ent_fine(r_ent,c_ent);
%     Cp(i)      = (Ent-H_old(i))/(T(i) - T_old(i));
    
    %Viscosity
    [~, r_vis] = min(abs(T(i)-t_vis_fine));
    [~, c_vis] = min(abs(P(i)-p_vis_fine));
    Mu(i)      = vis_fine(r_vis,c_vis);
end

end