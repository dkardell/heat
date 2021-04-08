function [D,G,I] = build_ops(Grid)

% author: Dominik A. Kardell
% modified from: Evan J. Ramos and Marc A. Hesse
% date: December 2020

% Description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.

% Input:
% Grid = structure containing all pertinent information about the grid.
% Grid.xmin
% Grid.xmax 
% Grid.Nx

% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);

% Output:
ex = ones(Grid.Nx,1);
Dx1 = spdiags([-ex ex]/Grid.dx, [0 1], Grid.Nx, Grid.Nx+1);
Ix = speye(Grid.Nx);
I = Ix;
if isfield(Grid,'ymin')
    ey = ones(Grid.Ny,1);
    Dy1 = spdiags([-ey ey]/Grid.dy, [0 1], Grid.Ny, Grid.Ny+1);
    Iy = speye(Grid.Ny);
    I = kron(Ix,Iy);
    Dx = kron(Dx1,Iy);
    Dy = kron(Ix,Dy1);
    D = horzcat(Dx,Dy);
    dof_f_bnd = [Grid.dof_f_xmin;Grid.dof_f_xmax;Grid.dof_f_ymin;Grid.dof_f_ymax];
    G = -D'; G(dof_f_bnd,:) = 0;
else
    D = Dx1;
    G = -D'; G(1,1) = 0; G(Grid.Nx+1,Grid.Nx) = 0;
end