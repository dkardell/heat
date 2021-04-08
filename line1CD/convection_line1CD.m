% File:    convection_line1AB.m
% Author:  Dominik A. Kardell, modified from Evan J. Ramos
% Date:    07 June 2020

clear, clc, close all

load phys_props.mat %Import physical property grids
dx = 37.5/(2^n_div); dy = dx; %[m] Grid spacing

%% Change grid size

dxn = 37.5; %New grid spacing in x-direction [m]
dyn = 37.5; %New grid spacing in y-direction [m]

xgrid           = 0:dx:(size(K_fine,2)-1)*dx;
ygrid           = 0:dy:(size(K_fine,1)-1)*dy;
[Xgrid,Ygrid]   = meshgrid(xgrid,ygrid);
xquery          = 0:dxn:max(xgrid);
yquery          = 0:dyn:max(ygrid);
[Xquery,Yquery] = meshgrid(xquery,yquery);

K_fine          = interp2(Xgrid,Ygrid,K_fine,Xquery,Yquery);
Ks_fine         = interp2(Xgrid,Ygrid,Ks_fine,Xquery,Yquery);
por_fine        = interp2(Xgrid,Ygrid,por_fine,Xquery,Yquery);
CPs_fine        = interp2(Xgrid,Ygrid,CPs_fine,Xquery,Yquery);
rho_fine        = interp2(Xgrid,Ygrid,rho_fine,Xquery,Yquery);
Vp_fine         = interp2(Xgrid,Ygrid,Vp_fine,Xquery,Yquery); sed_log = Vp_fine <= 3.5;

layer2A = (Vp_fine>=3 & Vp_fine<=5.16); 

% ocean_thick_vector        = interp1(xgrid,ocean_thick_vector,xquery)';
% ocean_thick_vector_smooth = interp1(xgrid,ocean_thick_vector_smooth,xquery)';

T_query = dyn/1000:dyn/1000:3.75;
T_x = [0,0.75,1.5,2.25,3,3.75];
T_values = 273.15+[2,32.38,64.4,96.28,128.14,160.09];
T_interp = interp1(T_x,T_values,T_query);

for i = 1:size(K_fine,2)
    for j = size(K_fine,1):-1:1
        if K_fine(j,i) > 1e-15
            crt(i) = j;
            T_bot(i) = T_interp(crt(i));
            break;
        end
    end
end

T_bot = T_interp(round(mean(crt)));
T_bot_vec = linspace(T_bot,T_bot,size(K_fine,2));
T_bot_vec = T_bot_vec';

dx = dxn; dy = dyn; %Update grid size

% figure(1)
% set(gcf,'units','normalized','outerposition',[0.5 0 0.5 0.9])
% subplot(4,1,1)
% contourf(por_fine,100,'LineStyle','none'); %axis equal;
% jetf = flipud(jet);
% set(gca, 'colormap', jetf); colorbar
% caxis([0 0.3])
% title('Porosity [frac]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,1,2)
% contourf(log10(K_fine),100,'LineStyle','none'); %axis equal;
% colormap(gca,jetf); colorbar 
% caxis([-20 -10])
% logvec = fliplr(logspace(-10,-20,6));
% colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
% title('Permeability [m^2]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,1,3)
% contourf(rho_fine,'LineStyle','none'); colorbar; %axis equal 
% xlabel('x (m)'); ylabel('z (m)');
% title('Density');
% 
% subplot(4,1,4)
% contourf(CPs_fine,'LineStyle','none'); colorbar; %axis equal 
% xlabel('x (m)'); ylabel('z (m)');
% title('Heat capacity');

%% Dimensional constants

%Space
Depth = size(K_fine,1)*dy;
Length = size(K_fine,2)*dx;

%Time
yrs2s    = 3600*24*365; %conversion factor to seconds from years
endtime  = 20000*yrs2s; %2000
t_totals = 5000; %200000 Total number of time steps saved
dt       = endtime/t_totals; %[s] time step
times    = (1:t_totals)*endtime/t_totals;

%% Build Grid and operators

Nx = size(K_fine,2); 
Ny = size(K_fine,1);

% Defining the computational grid and conductivities
Grid.xmin = 0; Grid.xmax = Length;  Grid.Nx = Nx;
Grid.ymin = 0; Grid.ymax = Depth;   Grid.Ny = Ny; 
Grid.Nz = 1;
Grid.psi_dir = 'xy';
Grid      = build_grid(Grid);
[Xc,Yc]   = meshgrid(Grid.xc,Grid.yc);
[Xf,Yf]   = meshgrid(Grid.xf,Grid.yf);
[Xfx,Yfx] = meshgrid(Grid.xf,Grid.yc);
[Xfy,Yfy] = meshgrid(Grid.xc,Grid.yf);

%Operators
[D,G,I]           = build_ops(Grid);

ocean_thick_vector = linspace(mean(ocean_thick_vector),mean(ocean_thick_vector),Nx)'; %3412,3318
ocean_thick_matrix = repmat(ocean_thick_vector',Ny,1);

%% Smooth data and plot

sy = 1; %Smoothing window size in x-direction (grid points)
sx = 0; %Smoothing window size in y-direction (grid points)

K_smooth   = smooth2a(K_fine,sy,sx);
Ks_smooth  = smooth2a(Ks_fine,sy,sx);
por_smooth = smooth2a(por_fine,sy,sx);
CPs_smooth = smooth2a(CPs_fine,sy,sx);
rho_smooth = smooth2a(rho_fine,sy,sx);

% figure(2)
% set(gcf,'units','normalized','outerposition',[0.3 0 0.7 0.9])
% subplot(4,2,1)
% contourf(por_fine, 100, 'LineStyle', 'none'); %axis equal;
% jetf = flipud(jet);
% set(gca, 'colormap', jetf); colorbar
% caxis([0 0.3])
% title('Porosity [frac]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,3)
% contourf(Xc,Yc,log10(K_fine),100,'LineStyle','none'); %axis equal;
% colormap(gca,jetf); colorbar 
% caxis([-20 -10])
% logvec = fliplr(logspace(-10,-20,6));
% colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
% title('Permeability [m^2]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,5)
% contourf(Xc,Yc,CPs_fine,'LineStyle','none'); colorbar; %axis equal 
% xlabel('x (m)'); ylabel('z (m)');
% title('Specific heat');
% 
% subplot(4,2,7)
% contourf(Xc,Yc,rho_fine,'LineStyle','none'); colorbar; %axis equal; 
% xlabel('x (m)'); ylabel('z (m)');
% title('Density');
% 
% subplot(4,2,2)
% contourf(por_smooth, 100, 'LineStyle', 'none'); %axis equal;
% set(gca, 'colormap', jetf); colorbar
% caxis([0 0.3])
% title('Smoothed Porosity [frac]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,4)
% contourf(Xc,Yc,log10(K_smooth),100,'LineStyle','none'); %axis equal;
% colormap(gca,jetf); colorbar 
% caxis([-20 -10])
% colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
% title('Smoothed Permeability [m^2]')
% xlabel('x [m]'); ylabel('z [m]');
% 
% subplot(4,2,6)
% contourf(Xc,Yc,CPs_smooth,'LineStyle','none'); colorbar; %axis equal; 
% xlabel('x (m)'); ylabel('z (m)');
% title('Smoothed specific heat');
% 
% subplot(4,2,8)
% contourf(Xc,Yc,rho_smooth,'LineStyle','none'); colorbar; %axis equal; 
% xlabel('x (m)'); ylabel('z (m)');
% title('Density');

K_fine   = K_smooth;
Ks_fine  = Ks_smooth;
por_fine = por_smooth;
CPs_fine = CPs_smooth;
rho_fine = rho_smooth;

%%

%Darcy
Mu       = 5e-5;      %[Pa*s] Dynamic viscosity of water at 20 deg C %%%%%%%%%%
rho_f    = 1030;        %[kg/m^3] Density of water
rho_s    = rho_fine;
grav     = 9.81;        %[m/s^2] Acceleration due to gravity on Earth
n        = 3;        %integer exponent for power law calc
tau      = sqrt(2);     %tortuosity of porosity medium, assuming packed spheres geometry
kvkh     = 1;

%Temperature
T_s       = 275.15;         %[K] surface temp
T_s_vec   = T_s*ones(Grid.Nx,1);

kf        = 0.65;        %[W/m/K] thermal conductivity of water
cp_f      = 4157;        %[J/(kg K)] specific heat of water

%% Create Medium

phi = por_fine; %por_smooth;
k = K_fine; %k_smooth;
ks = Ks_fine;

Kd = comp_mean(k,1,1,Grid);

cp_s = CPs_fine; %CPs_smooth;

%% Characteristic scaling terms

%Heat transport
kbar   = phi*kf + (1-phi).*ks; %phi/tau*kf + (1-phi).*ks;
rhobar = phi*rho_f*cp_f + (1-phi).*rho_s.*cp_s;

%% Boundary conditions pressure

% Parameters: pressure
Param.p.dof_dir   = Grid.dof_ymax;%[Grid.dof_ymin;Grid.dof_ymax]; %[Grid.dof_ymax];
Param.p.dof_f_dir = Grid.dof_f_ymax;%[Grid.dof_f_ymin;Grid.dof_f_ymax]; %[Grid.dof_f_ymax];
Param.p.g         = rho_f*grav*(ocean_thick_vector);%[rho_f*grav*(ocean_thick_vector+Depth);rho_f*grav*(ocean_thick_vector)]; %[rho_f*grav*(ocean_thick_vector)]; 
Param.p.dof_neu   = [];
Param.p.dof_f_neu = [];
Param.p.qb        = [];

[Bp,Np,fn_p] = build_bnd(Param.p,Grid,I);

%% Initialize pressure and temperature

%pressure
Ps                 = flipud(reshape(rho_f*grav*(ocean_thick_matrix+Yc),Grid.N,1));

%fluxes
%Q                = zeros(Grid.Nf,t_totals);

%Heat transport operators
theta_a          = 1;                          %advection solved explicitly
theta_d          = 0;                          %diffusion solved implicitly

Ts = zeros(Ny,Nx);
for i = 1:Nx
    Ts(:,i) = linspace(T_bot_vec(i),T_s,Ny);
end
Ts = reshape(Ts,[Nx*Ny 1]);

%temperature
T_i  = Ts;
T_ii = Ts;

%diagonalized matrices
Rhobar             = spdiags(rhobar(:),0,Grid.N,Grid.N);%%%%%%%%%%
Kbar               = comp_mean(kbar,1,kvkh,Grid); %%%%%%%%%%% 10 Aug 2017: arithmetic mean instead of harmonic mean
Rho_f              = rho_f;
KdMu               = Kd/Mu;

%Initialize loop variables
t = 0;
count = 1;

% mov = matfile('movie_1CD.mat','Writable',true);

%% Loop

while t < endtime
    
    t = t + dt;
    %% Solve steady state pressure --> Darcy Flux
    %pressure
    Lp             = -D*KdMu*G;
    fs_grav        = D*(KdMu*Rho_f*grav)*G*(ocean_thick_matrix(:)+Yc(:));
    fs_p           = zeros(Grid.N,1);
    fp             = fs_p + fn_p + fs_grav;
    P              = solve_lbvp(Lp,fp,Bp,Param.p.g,Np);
    
    q              = comp_flux_p(D,KdMu,Rho_f*grav*G*Yc(:),G,P,fs_p,Grid,Param.p);
    q(isnan(q)) = 0;
    
    qtop           = q(Grid.dof_f_ymax);
    qdown          = qtop(qtop<0);
    qdown_dof      = Grid.dof_ymax((qtop < 0));
    qdown_dof_f    = Grid.dof_f_ymax((qtop < 0));
        
    [qx,qy,qmags]  = center_flux(q,Grid);
    qy_up = ones(size(qy))*1e-20; qy_up(qy>1e-20) = qy(qy>1e-20);
    [psi,~,~]      = comp_streamfun(q,Grid);

    % 2nd order upwind of fluid fluxes
    Aq = flux_upwind(q,Grid);
    
    %% Boundary conditions temperature
    
    % Parameters: temperature
    Param.T.dof_dir   = [Grid.dof_ymin;qdown_dof]; %[Grid.dof_ymin;Grid.dof_ymax];
    Param.T.dof_f_dir = [Grid.dof_f_ymin;qdown_dof_f]; %[Grid.dof_f_ymin;Grid.dof_f_ymax];
    Param.T.g         = [T_bot_vec;T_s_vec((qtop<0))];

    Param.T.dof_neu   = [];%Grid.dof_ymin;
    Param.T.dof_f_neu = [];%Grid.dof_f_ymin;
    Param.T.qb        = [];%qt_b;

    [BT,NT,fn_T] = build_bnd(Param.T,Grid,I);
    
    %% Solve for temperature
    % Build heat operators
    Im_T             = @(theta_a,theta_d) ...
                     Rhobar + dt*(1-theta_a)*D*(rho_f*cp_f*Aq) ...
                            - dt*(1-theta_d)*D*(Kbar*G); %%%%%%
    Ex_T             = @(theta_a,theta_d) ...
                     Rhobar - dt*theta_a*D*(rho_f*cp_f*Aq)...
                            + dt*theta_d*D*(Kbar*G);     %%%%%%
    T_new      = solve_lbvp(Im_T(theta_a,theta_d),...
                            Ex_T(theta_a,theta_d)*T_i + ...
                            dt*(fn_T),BT,Param.T.g,NT);
   
    T_new(T_new<T_s) = T_s;
    T_new(T_new>T_bot) = T_bot;
                        
    %store for next time step
    T_ii = T_i;
    T_i  = T_new;
    
    %%
    T_ave(count) = mean(T_new-273.15);
    T_new_rect = reshape(T_new-273.15,Grid.Ny,Grid.Nx);
    T_top = T_new_rect(end,:);
    T_top_qup = T_top(qtop > 0); T_top_ave(count) = mean(T_top_qup); 
    sum_qup(count) = sum(qtop(qtop > 0)*dx);
    power_out = 1000 * cp_f * 1000 * (T_top-2).*qtop'; power_out(power_out<0) = 0; %milliWatts per m^2
    mean_power_out(count) = sum(1000 * cp_f * 1000 * ((T_top_qup-2).*qtop(qtop > 0)') / (length(qtop))); %milliWatts per m^2
    
    T_grad_top = (T_new_rect(end-1,:) - T_new_rect(end,:)) / dy;
    cond_power_out = 1000 * T_grad_top .* kbar(end,:); %[mW]
    % cond_power_out(cond_power_out<0) = 0;
    cond_mean_power_out(count) = sum(cond_power_out)/length(cond_power_out); %[mW]
    
    hfrat = power_out./(power_out + cond_power_out);
    hfrat(isnan(hfrat)) = 0;
        
    %% Update T-dependent variables
    %update fluid density
    [new_rho,~,new_mu] = steam_table(T_new,P,Grid);
    Rho_f         = comp_mean(reshape(new_rho,Grid.Ny,Grid.Nx),1,kvkh,Grid);
    new_Mu        = comp_mean(reshape(new_mu,Grid.Ny,Grid.Nx),1,kvkh,Grid);
    
    rho_plot = reshape(new_rho,Grid.Ny,Grid.Nx);
        
    KdMu = comp_mean(k./reshape(new_mu,Grid.Ny,Grid.Nx),1,kvkh,Grid);
    
    if xor(t == times(count),t > times(count))
        fprintf('\n%d of %d iterations\n',count,t_totals)
        count = count + 1;
    end
   
%% Plot while calculating
plotint = 100; %Plotting interval
if mod(count,plotint) == 1

figure(3)
    set(gcf,'units','normalized','outerposition',[0 0 0.7 1])
    subplot(3,1,1)
    contourf(Xc,Yc,T_new_rect,100,'LineStyle','none');
%     daspect([3 1 1]);
    colormap(jet);
    xlabel('Distance [m]')
    ylabel('Elevation [m]')
    c = colorbar;
    caxis([T_s-273.15 200]) %max(T_bot_vec)-273.15
    title(c,'T')
%     grid on
%     axis equal
    title(sprintf('Temperature after %5d years',round(times(count-1)/(yrs2s))))
%     set(gca,'xtick',[],'ytick',[]);

    subplot(3,1,2)
    contourf(Xc,Yc,log10(qmags),100,'LineStyle','none'); colorbar; hold on;
%     daspect([3 1 1]);
    logvec = fliplr(logspace(-6,-16,6));
    [qmax,qmaxd] = max(qmags); qmaxd = qmaxd*dy;
    jetf = flipud(jet);
    colormap(gca,jetf)
    colorbar('YTick', log10(logvec), 'YTickLabel', log10(logvec))
    caxis([-16 -6])
    contour(Xf,Yf,psi,20,'LineColor','w','LineWidth',0.5);
    plot(dx:dx:Length,qmaxd,'w.','MarkerSize',8); hold off;
    xlabel('x (m)'); ylabel('z (m)');
    title('Flow magnitude');
%     set(gca,'xtick',[],'ytick',[]);

    subplot(3,1,3)
    contourf(Xc,Yc,rho_plot,100,'LineStyle','none');
%     daspect([3 1 1]);
    colormap(gca,jetf);
    xlabel('Distance [m]')
    ylabel('Elevation [m]')
    c = colorbar;
%     caxis([T_s-273.15 200]) %max(T_bot_vec)-273.15
    title(sprintf('Density'))
%     set(gca,'xtick',[],'ytick',[]);
    
figure(4)
set(gcf,'units','normalized','outerposition',[0.7 0 0.3 1])
    subplot(5,1,1)
    plot(T_top)
    title('T along top boundary');
    subplot(5,1,2)
    bar(hfrat)
    xlim([0 length(hfrat)]); %ylim([0 0.6])
    title('Percentage adv. HF of total HF')
    subplot(5,1,3)
    plot(mean_power_out(2:end))
    title('Average advective power output [mW]');
    subplot(5,1,4)
    plot(cond_power_out,'r'); hold on
    plot(power_out,'b'); hold off
    ylim([0 1500])
    title('Heat flux along top bondary [mW]')
    subplot(5,1,5)
    plot(cond_mean_power_out(2:end))
    title('Average conductive power output [mW]');
    drawnow
   
end

%% Movie plot
% Xcp = Xc/1000;
% Ycp = (flipud(Yc)+mean(ocean_thick_vector))/1000;
% 
% plotint = 10; %Plotting interval
% if mod(count,plotint) == 1
%     figure(88)
%         set(gcf,'units','normalized','outerposition',[0.4 0.5 0.6 0.25])
%         contourf(Xcp,Ycp,T_new_rect,100,'LineStyle','none');
%         daspect([1 1 1]);
%         colormap(jet);
%         set(gca,'YDir','reverse')
%         xlabel('Distance [km]')
%         ylabel('Depth BSL [km]')
%         c = colorbar;
%         caxis([T_s-273.15 200]) %max(T_bot_vec)-273.15])
%         title(c,'^{o}C')
%         title(sprintf('31 Ma - Temperature after %5d years',round(times(count-1)/(yrs2s))))
%         drawnow
%         cdata = print('-RGBImage','-r85');
%         
%     F = im2frame(cdata);
%     % F = getframe(gcf);
%     mov = matfile('movie_1CD.mat','Writable',true);
%     mov.F(1,(count-1)/plotint) = F;
%  
% end

end

times = [0 times];

%% Calculate stats

% qup = qtop(qtop>0);
% 
% flux_sf = sum(qup)/length(qup);
% sprintf('Average flux across the seafloor = %d [m/s]',flux_sf)
% sprintf('FASF std dev = %d [m/s]',std(qup))
% sprintf('FASF std err = %d [m/s]',std(qup)/sqrt(length(qup)))
% % max_flow = max(qmax);
% sprintf('Average advective power output = %d mW/m^2',mean_power_out(end))
% sprintf('Advective std dev = %d mW/m^2',std(power_out))
% sprintf('Advective std err = %d mW/m^2',std(power_out)/sqrt(length(power_out)))
% sprintf('Average conductive power output = %d mW/m^2',cond_mean_power_out(end))
% sprintf('Conductive std dev = %d mW/m^2',std(cond_power_out))
% sprintf('Conductive std err = %d mW/m^2',std(cond_power_out)/sqrt(length(cond_power_out)))

% save('T_ave_line1CD.mat','T_ave');

%% Write out movie

% clear
% load movie_1CD.mat
% % F(1) = F(2);
% v = VideoWriter('convection_1CD.avi');
% open(v)
% writeVideo(v,F)
% close(v)