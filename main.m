% Clearing the workspace and adding all the necessary paths
clear; clc; close all

addpath(fullfile(cd,'Advection'));
addpath(fullfile(cd,'LagrangianDeform'));
addpath(fullfile(cd,'miscFtns'));

% % Load the data and generate a KDTree 
load(fullfile(cd,'Data','staticMesh.mat'));  Nt = size(time,2);
%% Parameters- NEEDS USER INPUT  

isStatic = 1; % 1- isStatic Mesh on the data, 0 - Meshdeforms/changes in every time
cpu_num = 2; % The number of parallel cpu cores to be alloted for advection 
% NOTE : The velocity data is copied across all cores therefore if velocity
% data is  y GB; The total RAM required for 'n' cores is n*y GB  

Nplot = 20;  % How many frames are being plotted on the videos generated in this code 
fntSz = 30; % FontSize of the plots 

% We compute the FTLE for trajectories from time(ct_i) -> time(ct_f)
% time array is loaded in along with the data 
ct_f = Nt; ct_i = 1; dt = 1*10^(-1);

%% Do some pre-compute for static Mesh and parallel pool Initiation

if isStatic == 1
    inputMesh = struct();
    
    xMesh = squeeze(x{1}); yMesh = squeeze(y{1}); zMesh = squeeze(z{1});
    inputMesh.faces = TrianT{1}; inputMesh.nodes = [squeeze(xMesh), squeeze(yMesh),squeeze(zMesh)];
    
    [face_mean_nodes,face_normals] = getFaceCenterAndNormals(inputMesh.faces,inputMesh.nodes);
    tree_model = KDTreeSearcher(face_mean_nodes);

    inputMesh.face_normals = face_normals; inputMesh.tree_model = tree_model;
    inputMesh.face_mean_nodes = face_mean_nodes;
end 

% Make parallel 
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)                                          % if parpool is not open
    parpool('local',cpu_num);
elseif (~isempty(poolobj)) && (poolobj.NumWorkers~=cpu_num)  % if parpool is not consistent with cpu_num
    delete(gcp)
    parpool('local',cpu_num);
end


%% See and visualize the data 

clc; close all; VideoFileName=fullfile(cd,'saveResults','vizVelocity'); coarse = 2; % How much to coarsen velocity field
writerObj = VideoWriter(VideoFileName,'MPEG-4');writerObj.FrameRate = 5;  writerObj.Quality = 100;
writerObj.FileFormat; open(writerObj); count = 0;

figure('Color','w','units','normalized','outerposition',[0 0 1 1]);

for idXMesh=1:round(Nt/Nplot):Nt
    count = count+1;
    t_i = time(idXMesh);

    v1Data = squeeze(v{1,idXMesh}); v2Data = squeeze(v{2,idXMesh}); v3Data = squeeze(v{3,idXMesh});
    xData = squeeze(x{idXMesh}); yData = squeeze(y{idXMesh}); zData = squeeze(z{idXMesh});
    TriData = TrianT{idXMesh}; vMag_Data = sqrt(v1Data.^2+v2Data.^2+v3Data.^2);

    % Plot the quantities
    subplot(1,1,1)
    trisurf(TriData,xData,yData,zData,vMag_Data,'FaceAlpha',1,'Edgecolor','none')
    colormap(parula); shading interp; hold on 

    quiver3(xData(1:coarse:end),yData(1:coarse:end),zData(1:coarse:end),...
        v1Data(1:coarse:end),...
        v2Data(1:coarse:end),...
        v3Data(1:coarse:end),...
        'LineWidth',1,'Color','k','MaxHeadSize',20); colorbar;

    sgtitle(sprintf('t = %.2f',t_i),'Interpreter','latex','FontSize',fntSz);
    daspect([1 1 1]); axis tight off;  hold off
    
    set(gca,'FontSize',fntSz)
    mov(count) = getframe(gcf);       
end 

writeVideo(writerObj,mov); close(writerObj); 

%% Advect the grid backward 
clc; close;

tau0_Advect = 0; tauf_Advect = time(ct_f)-time(ct_i); 
xqBck = x{ct_f}; yqBck = y{ct_f}; zqBck = z{ct_f}; % Initial condition for backward advect 
xqFor = x{ct_i}; yqFor = y{ct_i}; zqFor = z{ct_i}; % Initial condition for forward advect

tf_Absolute = time(ct_f);
tauArr_Advect = tau0_Advect:dt:tauf_Advect;
tauSave_Advect = tauArr_Advect(1:1:end); 

if isStatic == 1
    [xt_Advect,yt_Advect,zt_Advect] = AdvctBckStatic(cpu_num,tf_Absolute,tauSave_Advect, ...
        tauArr_Advect,xqBck,yqBck,zqBck,time,v,inputMesh);
else
    [xt_Advect,yt_Advect,zt_Advect] = AdvctBck(cpu_num,tf_Absolute,tauSave_Advect, ...
        tauArr_Advect,xqBck,yqBck,zqBck,time,v,x,y,z,TrianT);
end 

% % Advect the grid forward
t0_Advect_for = time(ct_i); tf_Advect_for = time(ct_f);
tArr_Advect_for = t0_Advect_for:dt:tf_Advect_for;
tSave_Advect_for = tArr_Advect_for(1:1:end); 

%Can't save all the data, will save at selective times
if isStatic == 1
    [xt_Advect_for,yt_Advect_for,zt_Advect_for] = AdvctForStatic(cpu_num,tSave_Advect_for, ...
        tArr_Advect_for,xqFor,yqFor,zqFor,time,v,inputMesh);
else
    [xt_Advect_for,yt_Advect_for,zt_Advect_for] = AdvctFor(cpu_num,tSave_Advect_for, ...
        tArr_Advect_for,xqFor,yqFor,zqFor,time,v,x,y,z,TrianT);
end 

%% Plot the advection

% % Backward advection
clc; close all; VideoFileName=fullfile(cd,'saveResults','bckAdvct'); coarse = 10; % How much to coarsen velocity field
writerObj = VideoWriter(VideoFileName,'MPEG-4');writerObj.FrameRate = 5;  writerObj.Quality = 100;
writerObj.FileFormat; open(writerObj); count = 0;

figure('Color','w','units','normalized','outerposition',[0 0 1 1])

for i=1:round(numel(tauSave_Advect)/Nplot):numel(tauSave_Advect)
    count = count+1;

    t_i = tf_Absolute-tauSave_Advect(i);

    % Load the mesh quantities 
    [~,idXMesh] = min(abs(time-t_i));
    v1Data = -squeeze(v{1,idXMesh}); v2Data = -squeeze(v{2,idXMesh}); v3Data = -squeeze(v{3,idXMesh});
    xData = squeeze(x{idXMesh}); yData = squeeze(y{idXMesh}); zData = squeeze(z{idXMesh});
    TriData = TrianT{idXMesh}; vMag_Data = sqrt(v1Data.^2+v2Data.^2+v3Data.^2);
    %Load the location of the advected particles 
    xPts_i = xt_Advect(i,:); yPts_i = yt_Advect(i,:); zPts_i = zt_Advect(i,:);

    % Plot the quantities
    trisurf(TriData,xData,yData,zData,vMag_Data,'FaceAlpha',1,'Edgecolor','none')
    colormap(parula); shading interp; hold on 
    quiver3(xData(1:coarse:end),yData(1:coarse:end),zData(1:coarse:end),...
        v1Data(1:coarse:end),...
        v2Data(1:coarse:end),...
        v3Data(1:coarse:end),...
        'LineWidth',1,'Color','k','MaxHeadSize',20)
    scatter3(xPts_i,yPts_i,zPts_i,10,'r','filled'); hold off 
    title(sprintf('t = %.2f',t_i),'Interpreter','latex');daspect([1 1 1]); axis tight off
    
    set(gca,'FontSize',fntSz)
    mov(count) = getframe(gcf);       
end 
writeVideo(writerObj,mov); close(writerObj); 

% % Forward Advection Video
clc; close all; VideoFileName=fullfile(cd,'saveResults','forAdvct'); coarse = 10; % How much to coarsen velocity field
writerObj = VideoWriter(VideoFileName,'MPEG-4');writerObj.FrameRate = 5;  writerObj.Quality = 100;
writerObj.FileFormat; open(writerObj); count = 0;
figure('Color','w','units','normalized','outerposition',[0 0 1 1])

for i=1:round(numel(tSave_Advect_for)/Nplot):numel(tSave_Advect_for)
    count = count+1;

    t_i = tSave_Advect_for(i);

    % Load the mesh quantities 
    [~,idXMesh] = min(abs(time-t_i));
    v1Data = squeeze(v{1,idXMesh}); v2Data = squeeze(v{2,idXMesh}); v3Data = squeeze(v{3,idXMesh});
    xData = squeeze(x{idXMesh}); yData = squeeze(y{idXMesh}); zData = squeeze(z{idXMesh});
    TriData = TrianT{idXMesh}; vMag_Data = sqrt(v1Data.^2+v2Data.^2+v3Data.^2);
    %Load the location of the advected particles 
    xPts_i = xt_Advect_for(i,:); yPts_i = yt_Advect_for(i,:); zPts_i = zt_Advect_for(i,:);

    % Plot the quantities
    trisurf(TriData,xData,yData,zData,vMag_Data,'FaceAlpha',1,'Edgecolor','none')
    colormap(parula); shading interp; hold on 
    quiver3(xData(1:coarse:end),yData(1:coarse:end),zData(1:coarse:end),...
        v1Data(1:coarse:end),...
        v2Data(1:coarse:end),...
        v3Data(1:coarse:end),...
        'LineWidth',1,'Color','k','MaxHeadSize',20)
    scatter3(xPts_i,yPts_i,zPts_i,10,'r','filled'); hold off 
    title(sprintf('t = %.2f',t_i),'Interpreter','latex'); daspect([1 1 1]); axis tight off
    
    set(gca,'FontSize',fntSz)
    mov(count) = getframe(gcf);       
end 
writeVideo(writerObj,mov); close(writerObj); 

%% Calculate and Visualize the FTLE values 
clc; close all; viewAngle = [-107,20];

% Compute bw-FTLE

x0 = squeeze(xt_Advect(1,:)); y0 = squeeze(yt_Advect(1,:)); z0 = squeeze(zt_Advect(1,:));
xf = squeeze(xt_Advect(end,:)); yf = squeeze(yt_Advect(end,:)); zf = squeeze(zt_Advect(end,:));
tf = time(ct_i); ti = time(ct_f); % Lagrangian time interval of analysis 

Tri_pts_Uni = TrianT{ct_f}; % The triangulation of the mesh at tf
xfData = squeeze(x{1}); yfData = squeeze(y{1}); zfData = squeeze(z{1}); Trif_Data = squeeze(TrianT{1});
[lambdaField,lambdaIsoField,maxDefEigVec] = lagDefCompute(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_pts_Uni,Trif_Data,tf,ti);

f= figure('color','w','Units','normalized','OuterPosition',[0.0060 0.2435 0.9878 0.619]);

subplot(1,4,1)
trisurf(Tri_pts_Uni,x0,y0,z0,lambdaIsoField,'FaceAlpha',1,'Edgecolor','none');
hold on;  colorbar('southoutside'); axis equal; shading interp;  axis off;
title(sprintf('$ {}_{iso}\\Lambda_{%.2f}^{%.2f} $',time(ct_f),time(ct_i)),'Interpreter','latex'); set(gca,'FontSize',fntSz);
view(viewAngle); camva(7);

subplot(1,4,2)
trisurf(Tri_pts_Uni,x0,y0,z0,lambdaField,'FaceAlpha',1,'Edgecolor','none'); hold on 
quiver3(x0,y0,z0,maxDefEigVec(1,:),maxDefEigVec(2,:),maxDefEigVec(3,:),'k','Linewidth',1,'ShowArrowHead','off'); hold off
axis equal; shading interp; axis off; view(viewAngle); camva(7);
title(sprintf('$ \\Lambda_{%.2f}^{%.2f},\\mathbf{\\zeta}_{%.2f}^{%.2f} $',time(ct_f),time(ct_i),time(ct_f),time(ct_i)),'Interpreter','latex'); set(gca,'FontSize',fntSz);
colorbar('southoutside')

% Compute FW-FTLE
x0 = squeeze(xt_Advect_for(1,:)); y0 = squeeze(yt_Advect_for(1,:)); z0 = squeeze(zt_Advect_for(1,:));
xf = squeeze(xt_Advect_for(end,:)); yf = squeeze(yt_Advect_for(end,:)); zf = squeeze(zt_Advect_for(end,:));
tf = time(ct_f); ti = time(ct_i); % Lagrangian time interval of analysis 

Tri_pts_Uni = TrianT{1}; xfData = squeeze(x{end}); yfData = squeeze(y{end}); zfData = squeeze(z{end}); Trif_Data = squeeze(TrianT{end});

[lambdaField,lambdaIsoField,maxDefEigVec] = lagDefCompute(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_pts_Uni,Trif_Data,tf,ti);

subplot(1,4,3)
trisurf(Tri_pts_Uni,x0,y0,z0,lambdaIsoField,'FaceAlpha',1,'Edgecolor','none');
colorbar('southoutside'); axis equal; shading interp;  axis off;
view(viewAngle); camva(7);
title(sprintf('$ {}_{iso}\\Lambda_{%.2f}^{%.2f} $',time(ct_i),time(ct_f)),'Interpreter','latex'); set(gca,'FontSize',fntSz);

subplot(1,4,4)
trisurf(Tri_pts_Uni,x0,y0,z0,lambdaField,'FaceAlpha',1,'Edgecolor','none');hold on 
quiver3(x0,y0,z0,maxDefEigVec(1,:),maxDefEigVec(2,:),maxDefEigVec(3,:),'k','Linewidth',1,'ShowArrowHead','off'); hold off
colorbar('southoutside'); axis equal; shading interp;  axis off;
view(viewAngle); camva(7);
title(sprintf('$ \\Lambda_{%.2f}^{%.2f},\\mathbf{\\zeta}_{%.2f}^{%.2f} $',time(ct_i),time(ct_f),time(ct_i),time(ct_f)),'Interpreter','latex'); set(gca,'FontSize',fntSz);


