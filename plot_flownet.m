function [] = plot_flownet(Nh,Ns,h,PSI,h_style,psi_style,Grid)

% author: Dominik A. Kardell
% modified from: Evan J. Ramos and Marc A. Hesse 
% date: December 2020

% Description: Plots a flownet with an Nh equally spaced head contours
% and Ns equally spaced streamlines.

% Input: Nh = number of head contours
% Ns = number of streamlines
% h = Grid.N by 1 column vector of heads
% PSI = Ny+1 by Nx+1 matrix containing the stream function
% h_style = string specifying the linestyle for the head contours
% psi_style = string specifying the linestyle for the streamlines
% Grid = structure containing all information about the grid.

% Output: 
x = linspace(0,Grid.xmax,Grid.Nx); z = linspace(0,Grid.ymax,Grid.Ny);
[X,Z] = meshgrid(x,z);
a = linspace(0,Grid.xmax,Grid.Nx+1); b = linspace(0,Grid.ymax,Grid.Ny+1);
[A,B] = meshgrid(a,b);

% figure(fig)
% contour(X,Z,hn,Nh,h_style)
% hold on
contour(A,B,PSI,Ns,psi_style)
axis equal
xlabel('x (m)')
ylabel('z (m)')