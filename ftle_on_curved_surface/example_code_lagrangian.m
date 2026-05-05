% This is an example of how to use the code 
clear; clc; close all;

% Loading dataset and functions
load('./example_data.mat')
addpath('./functions/')
Nt = numel(mesh_time);

% Parameters - need user input
delta = 0.2; % geodesic distance over which to compute deformation
ct0 = 1; ctf = Nt; % time-interval of analysis 

% Visualize the velocity field 
for ct = 1:Nt
    xct = mesh_r{ct,1};
    yct = mesh_r{ct,2};
    zct = mesh_r{ct,3}; 

    v1ct = mesh_v{ct,1};
    v2ct = mesh_v{ct,2};
    v3ct = mesh_v{ct,3};

    vMag = sqrt(v1ct.^2+v2ct.^2+v3ct.^2);

    trisurf(mesh_F{ct},xct,yct,zct,vMag,'Edgecolor','none'); colorbar
    hold on; quiver3(xct,yct,zct,v1ct,v2ct,v3ct,'k'); hold off
    axis equal off
    title(sprintf('v(x, t= %.1f) ',mesh_time(ct))...
        ,'FontSize',24,'Interpreter','latex')
    shading interp
    pause(0.2)
end 
% Splicing the dataset to corresponding dataset
% Performs analysis from mesh_time(ct0) to mesh_time(ctf)
mesh_F = mesh_F(ct0:ctf,1);
mesh_time = mesh_time(ct0:ctf);
mesh_v = mesh_v(ct0:ctf,:);
mesh_r = mesh_r(ct0:ctf,:);

% Advect the particles forward 
rt_arr = advect_for(mesh_F,mesh_time,mesh_v,mesh_r);

% Compute the deformation
mesh_r0 = cell2mat(mesh_r(1,:)); % Mesh at 
mesh_rf = cell2mat(mesh_r(end,:)); % Final particle positions after advection
mesh_F0 = cell2mat(mesh_F(1)); % Initial mesh faces
mesh_Ff = cell2mat(mesh_F(end)); % Final mesh faces
rf = cell2mat(rt_arr(end,:));

% Compute the deformation metric
[L2,L1,V0,Vf] = compute_deform(mesh_r0,mesh_rf,mesh_F0,mesh_Ff,rf,delta);
FTLE = 0.5*log(L2)./(mesh_time(end)-mesh_time(1));
J = 0.5*(log(L1)+log(L2))./(mesh_time(end)-mesh_time(1));
aniso = 0.5*(log(L2)-log(L1));

%% Visualize the deformation metric
close all; fontSz = 24; camAmp = 8; cam_view = [-33 1.7];

% Define the custom dark background color
bgColor = [13 13 13]/255; 

% Set the figure background color
f1 = figure('Position',[72 767 1609 426], 'Color', bgColor);

% --- Subplot 1 ---
ax1 = subplot(1,3,1);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),J,'Edgecolor','none')
c1 = colorbar; 
c1.Color = 'w'; % Make colorbar text white
axis equal off; 
ax1.FontSize = fontSz; 
ax1.Color = bgColor;          % Axes background color
ax1.XColor = 'w'; ax1.YColor = 'w'; ax1.ZColor = 'w'; % Axes text/ticks color
camva(camAmp); shading interp; 
t1 = title(sprintf('$$ {}_{iso}\\Lambda_{%.1f}^{%.1f }(\\mathbf{x}_0)$$',mesh_time(1),mesh_time(end)),...
    'FontSize',fontSz,'Interpreter','latex'); 
t1.Color = 'w'; % Make title text white
view(cam_view)

% --- Subplot 2 ---
ax2 = subplot(1,3,2);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),aniso,'Edgecolor','none')
c2 = colorbar; 
c2.Color = 'w'; 
axis equal off;
ax2.FontSize = fontSz; 
ax2.Color = bgColor;
ax2.XColor = 'w'; ax2.YColor = 'w'; ax2.ZColor = 'w';
camva(camAmp); shading interp; 
t2 = title(sprintf('$${}_{aniso}\\Lambda_{%.1f}^{%.1f }(\\mathbf{x}_0)$$',mesh_time(1),mesh_time(end)),...
    'FontSize',fontSz,'Interpreter','latex'); 
t2.Color = 'w';
view(cam_view)

% --- Subplot 3 ---
ax3 = subplot(1,3,3);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),FTLE,'Edgecolor','none')
idx = datasample(1:size(V0,1),1000,'Replace',false); % Undersampling the eigenvectors to display
hold on; 
% Changed quiver color from 'k' to 'w' for visibility on dark background
quiver3(mesh_r0(idx,1),mesh_r0(idx,2),mesh_r0(idx,3),...
    V0(idx,1),V0(idx,2),V0(idx,3),'k','ShowArrowHead','off'); 
hold off
c3 = colorbar; 
c3.Color = 'w'; 
axis equal off;
ax3.FontSize = fontSz; 
ax3.Color = bgColor;
ax3.XColor = 'w'; ax3.YColor = 'w'; ax3.ZColor = 'w';
camva(camAmp); shading interp; 
t3 = title(sprintf('$ \\Lambda_{%.1f}^{%.1f }(\\mathbf{x}_0), \\xi_{%.1f}^{%.1f }(\\mathbf{x}_0)$', ...
    mesh_time(1), mesh_time(end), mesh_time(1), mesh_time(end)), ...
    'FontSize', fontSz, 'Interpreter', 'latex');
t3.Color = 'w';
view(cam_view)

