% Example code demonstrating the pipeline to compute the Lagrangian
% deformation for 2D/3D Euclidean flows 

clc; clear; close all

dxy = 0.01; 
x = 0:dxy:2; y = 0:2*dxy:1;
[X,Y] = meshgrid(x,y); [N1,N2] = size(X);

% % Advect tracer particles 
r0 = [X(:);Y(:)];
options = odeset('AbsTol', 1e-8, 'RelTol', 1e-6);
[t,r_t] = ode45(@(t,y)velode45(t,y),0:0.1:1,r0,options);
num_steps = length(t); num_ICs = size(r0,1)/2;
r_t = reshape(r_t.',N1,N2,2, num_steps);

%% Compute Lagrangian deformation  

x0 = squeeze(r_t(:,:,1,1));
y0 = squeeze(r_t(:,:,2,1));

xf = squeeze(r_t(:,:,1,end));
yf = squeeze(r_t(:,:,2,end));

delta = 0.03; % The spatial averaging parameter 
r0 = [x0(:),y0(:)]; rf = [xf(:),yf(:)];
[sv_list,vec0_list] = compute_deform_euclidean(r0,rf,delta);

%% Visualize the result
figure('Position',[1583 536 675 547])
fntSz = 24; cr = 2;

% Reshape the velocity field 
FTLE = log(sv_list(:,end))/(max(t)-min(t));
FTLE = reshape(FTLE,size(X));
eigvec = squeeze(vec0_list(:,:,end));

contourf(X,Y,FTLE,'EdgeColor','none'); c = colorbar; 
hold on; quiver(r0(1:cr:end,1),r0(1:cr:end,2),...
    eigvec(1:cr:end,1),eigvec(1:cr:end,2),...
    1,'r',"ShowArrowHead","off"); hold off
xlim([0,2]); ylim([0,1]); axis equal
ax = gca; ax.FontSize = fntSz; ax.LineWidth = 2;
c.FontSize = fntSz; c.LineWidth = 2;
title(sprintf('$$ \\Lambda_{%.1f}^{%.1f},\\xi_{%.1f}^{%.1f} $$',...
    min(t),max(t),min(t),max(t)),'Interpreter','latex','FontSize',fntSz)

%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%
function dydt = velode45(t,y)
dydt = 0*y;

Nq = size(y,1)/2;
xq = y(1:Nq); yq = y(Nq+1:end);
[u,v] = velocity(xq,yq,t);

dydt(1:Nq) = u;
dydt(Nq+1:end) = v;

end 

function psi = stream_func(x,y,t)
g_t = 0.3*sin(4*pi*t)+0.1*sin(2*t);
psi = sin(pi*(x-g_t)).*sin(pi*y);
end


function [u,v] = velocity(X,Y,t)
% The streamfunction for a 2D model of a Rayleigh bernard convection cell
dxy = 10^(-4);

dPsi_dx = (stream_func(X+dxy,Y,t)-stream_func(X,Y,t))./dxy;
dPsi_dy = (stream_func(X,Y+dxy,t)-stream_func(X,Y,t))./dxy;

u =  dPsi_dy;  % u = dPsi/dy
v = -dPsi_dx;  % v = -dPsi/dx

end 
