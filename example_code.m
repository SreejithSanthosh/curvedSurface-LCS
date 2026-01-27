% This is an example of how to use the code 
clear; clc; close all;

% Loading dataset and functions
load('./example_data.mat')
addpath('./functions/')

% Visualize the velocity field 
Nt = numel(mesh_time);

for ct = 1:Nt
    xct = mesh_r{ct,1};
    yct = mesh_r{ct,2};
    zct = mesh_r{ct,3}; 

    v1ct = mesh_v{ct,1};
    v2ct = mesh_v{ct,2};
    v3ct = mesh_v{ct,3};
    
    vMag = sqrt(v1ct.^2+v2ct.^2+v3ct.^2);

    trisurf(mesh_F{ct},xct,yct,zct,vMag)
    hold on; quiver3(xct,yct,zct,v1ct,v2ct,v3ct); hold off
    axis equal off
end 

% Advect the particles forward 
rt_arr = advect_for(mesh_F,mesh_time,mesh_v,mesh_r);

% Compute the deformation

% Define the mesh from 0->ct 
ct = numel(mesh_time); % Select the index to compute deformation 


mesh_r0 = cell2mat(mesh_r(1,:)); % Mesh at 
mesh_rf = cell2mat(mesh_r(ct,:)); % Final particle positions after advection
mesh_F0 = cell2mat(mesh_F(1)); % Initial mesh faces
mesh_Ff = cell2mat(mesh_F(ct)); % Final mesh faces
rf = cell2mat(rt_arr(end,:));

% Compute the deformation metric
delta = 0.2;
[L2,L1,V0,Vf] = compute_deform(mesh_r0,mesh_rf,mesh_F0,mesh_Ff,rf,delta);

FTLE = L2./(mesh_time(ct)-mesh_time(1));
J = 0.5*(L1+L2)./(mesh_time(ct)-mesh_time(1));

%% Visualize the deformation metric
close all; fontSz = 24; camAmp = 8; cam_view = [-33 1.7];
f1 = figure('Position',[102 192 638 753]);
theme('light')
ax = subplot(2,1,1);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),J,'Edgecolor','none')
colorbar; axis equal off; ax.FontSize = fontSz; camva(camAmp); 
shading interp; title(sprintf('$$ {}_{iso}\\Lambda_{0}^{%.1f }(\\mathbf{x}_0)$$',mesh_time(ct)),...
    'FontSize',fontSz,'Interpreter','latex'); view(cam_view)

ax = subplot(2,1,2);
trisurf(mesh_F0,mesh_r0(:,1),mesh_r0(:,2),mesh_r0(:,3),FTLE,'Edgecolor','none')
idx = datasample(1:size(V,1),1000,'Replace',false); % Undersampling the eigenvectors to display
hold on; quiver3(mesh_r0(idx,1),mesh_r0(idx,2),mesh_r0(idx,3),...
    V(idx,1),V(idx,2),V(idx,3),'k','ShowArrowHead','off'); hold off
colorbar; axis equal off;ax.FontSize = fontSz; camva(camAmp); 
shading interp; title(sprintf('$$\\Lambda_{0}^{%.1f }(\\mathbf{x}_0)$$',mesh_time(ct)),...
    'FontSize',fontSz,'Interpreter','latex'); view(cam_view)