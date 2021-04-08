function [Grid] = build_grid(Grid)

%% 

% author: Dominik A. Kardell
% modified from: Evan J. Ramos and Marc A. Hesse
% date: December 2020

% Description: This function computes takes in minimal definition of the 
% computational domain and grid and computes all containing all pertinent 
% information about the grid. 

% Input: 
% Grid.xmin = x-min boundary of the domain 
% Grid.xmax = x-max boundary of the domain 
% Grid.Nx = number of grid cells: x-direction

% Grid.ymin = y-min boundary of the domain
% Grid.ymax = y-max boundary of the domain
% Grid.Ny = number of grid cells: y-direction

% Grid.dz = z-width that the user can input

% Output:
% In this case, I'll assign dummy variables to the given fields if there is
% only a 1D case. ***Arbitrary values will be assigned***

if ~isfield(Grid,'ymin')
    Grid.ymin = 0;
    Grid.ymax = 10;
    Grid.Ny = 1;
    Grid.D2 = false; % Field that exists only if Grid is for 1D flow
                     % **** to be used in build_ops.m function ****
end

Grid.Lx = Grid.xmax - Grid.xmin; % length of the domain
Grid.Ly = Grid.ymax - Grid.ymin; % depth of the domain

Grid.dx = Grid.Lx/Grid.Nx; % cell width
Grid.dy = Grid.Ly/Grid.Ny; % cell length

mid_cell_x = Grid.dx/2;
Grid.xc = ((Grid.xmin+mid_cell_x):Grid.dx:(Grid.xmax-mid_cell_x))'; % vector of cell center locations: x-direction
if Grid.xmin < 0 %assuming it is symmetric about zero
    vec = mid_cell_x:Grid.dx:Grid.xmax-mid_cell_x;
    Grid.xc = [-fliplr(vec) vec]';
end
mid_cell_y = Grid.dy/2;
Grid.yc = ((Grid.ymin+mid_cell_y):Grid.dy:(Grid.ymax-mid_cell_y))'; % vector of cell center locations: y-direction

Grid.xf = (Grid.xmin:Grid.dx:Grid.xmax)'; % vector of cell face locations: x-direction
Grid.yf = (Grid.ymin:Grid.dy:Grid.ymax)'; % vector of cell face locations: y-direction

Grid.N = Grid.Nx*Grid.Ny; % total number of cells
Grid.Nfx = (Grid.Nx+1)*Grid.Ny; % number of fluxes in x-direction 
Grid.Nfy = (Grid.Ny+1)*Grid.Nx; % number of fluxes in y-direction 
Grid.Nf = Grid.Nfx + Grid.Nfy; % total number of fluxes

if xor(Grid.Nx == 1,Grid.Ny == 1) %1D
    Grid.Nf = Grid.Nfx;
end

%% Degrees of Freedom

% Grid.dof = vector from 1 to N containing the degrees of freedom
%            i.e. cell numbers 
Grid.dof = 1:Grid.N;

% Grid.dof_f = vector from 1 to Nf containing the degrees of freedom
%              associated with cell faces
Grid.dof_f = 1:Grid.Nf;

% Creating various matrices to ease in the selection of the degrees of
% freedom pertaining to cell centers and cell faces;

cell_center_mat = flipud(reshape(Grid.dof,Grid.Ny,Grid.Nx));
x_flux_mat = flipud(reshape(1:Grid.Nfx,Grid.Ny,Grid.Nx+1));

Grid.dof_xmin = flipud(cell_center_mat(:,1)); % degrees of freedom corresponding to the left boundary
Grid.dof_xmax = flipud(cell_center_mat(:,end)); % degrees of freedom corresponding to the right boundary 
Grid.dof_ymin = cell_center_mat(end,:)'; % degrees of freedom corresponding to the bottom boundary
Grid.dof_ymax = cell_center_mat(1,:)'; % degrees of freedom corresponding to the top boundary

if xor(Grid.Nx == 1,Grid.Ny == 1) %1D
    Grid.dof_xmin = 1;
    Grid.dof_xmax = Grid.Nx;
end

Grid.dof_f_xmin = flipud(x_flux_mat(:,1)); % The face/flux number on the left boundary
Grid.dof_f_xmax = flipud(x_flux_mat(:,end)); % The face/flux number on the right boundary
Grid.dof_f_ymin = Grid.Nfx + (1:Grid.Ny+1:Grid.Nfy-Grid.Ny+1)'; % The face/flux number on the bottom boundary
Grid.dof_f_ymax = Grid.Nfx + (Grid.Ny+1:Grid.Ny+1:Grid.Nfy)'; % The face/flux number on the top boundary

if xor(Grid.Nx == 1,Grid.Ny == 1) %1D
    Grid.dof_f_xmin = 1;
    Grid.dof_f_xmax = Grid.Nfx;
end

%% Volume elements

% Establishing arbitrary z-coordinate values (widths)
if ~isfield(Grid,'dz')
    Grid.dz = 1;
end

Grid.V = (Grid.dx*Grid.dy*Grid.dz)*ones(Grid.Nf,1); % A column vector of length Nf that contains the cell volumes
Grid.A = ones(Grid.Nf,1); % A column vector of length Nf that containts all face areas
% boundaries
Grid.A(1:Grid.Nfx) = Grid.dy*Grid.dz;
Grid.A(Grid.Nfx+(1:Grid.Nfy)) = Grid.dx*Grid.dz;

if (Grid.Nx == 1 && Grid.Ny > 1) || (Grid.Nx > 1 && Grid.Ny == 1) %1D
    Grid.A = Grid.dx*Grid.dz*ones(Grid.Nf,1);
    Grid.V = Grid.dx^2*Grid.dz*ones(Grid.Nf,1);
end


%% Stream function elements

if ~isfield(Grid,'psi_x0')
    Grid.psi_x0 = 'xmin_ymin';
end

if ~isfield(Grid,'psi_dir')
    Grid.psi_dir = 'yx';
end

end