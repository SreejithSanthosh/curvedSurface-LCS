% We are wrting a more refined version of the code 
clear; clc; close all

load('./formatData/dataHeart.mat')
Nt = size(x,1); isStatic = 0; % 1- isStatic Mesh on the data, 0 - Meshdeforms/changes in every time
Nplot = 20;  % How many frames are being plotted on the videos 

% Define all the properties of the mesh if Static to 
if isStatic == 1
    inputMesh = struct();
    
    xMesh = squeeze(x{1}); yMesh = squeeze(y{1}); zMesh = squeeze(z{1});
    inputMesh.faces = TrianT{1}; inputMesh.nodes = [squeeze(xMesh), squeeze(yMesh),squeeze(zMesh)];
    
    [face_mean_nodes,face_normals] = getFaceCenterAndNormals(inputMesh.faces,inputMesh.nodes);
    tree_model = KDTreeSearcher(face_mean_nodes);

    inputMesh.face_normals = face_normals; inputMesh.tree_model = tree_model;
    inputMesh.face_mean_nodes = face_mean_nodes;
end 

%% See and visualize the data 

clc; close all; VideoFileName=strcat(cd,sprintf('/vizVelocity')); coarse = 10; % How much to coarsen velocity field
writerObj = VideoWriter(VideoFileName,'MPEG-4');writerObj.FrameRate = 5;  writerObj.Quality = 100;
writerObj.FileFormat; open(writerObj); count = 0;
figure('Color','w','units','normalized','outerposition',[0 0 1 1])

for idXMesh=1:round(Nt/Nplot):Nt
    count = count+1;
    t_i = time(idXMesh);

    v1Data = squeeze(v{1,idXMesh}); v2Data = squeeze(v{2,idXMesh}); v3Data = squeeze(v{3,idXMesh});
    xData = squeeze(x{idXMesh}); yData = squeeze(y{idXMesh}); zData = squeeze(z{idXMesh});
    TriData = TrianT{idXMesh}; vMag_Data = sqrt(v1Data.^2+v2Data.^2+v3Data.^2);

    % Plot the quantities
    trisurf(TriData,xData,yData,zData,vMag_Data,'FaceAlpha',1,'Edgecolor','none')
    colormap(parula); shading interp; hold on 
    quiver3(xData(1:coarse:end),yData(1:coarse:end),zData(1:coarse:end),...
        v1Data(1:coarse:end),...
        v2Data(1:coarse:end),...
        v3Data(1:coarse:end),...
        'LineWidth',1,'Color','k','MaxHeadSize',20)
    title(sprintf('t = %.2f',t_i));daspect([1 1 1]); axis tight;  hold off

    mov(count) = getframe(gcf);       
end 

writeVideo(writerObj,mov); close(writerObj); 


%% Construct parallel pool.
% EDIT: Set number cores.
MaxNcores = 20; NpDG = 20;
cpu_num = MaxNcores;

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)                                          % if parpool is not open
    parpool('local',cpu_num);
elseif (~isempty(poolobj)) && (poolobj.NumWorkers~=cpu_num)  % if parpool is not consistent with cpu_num
    delete(gcp)
    parpool('local',cpu_num);
end

%% Specify a loop to go over bFTLE  
clc; close; ct_fArr = 2:1:Nt; % This is the time points you store the velocity values 
time = 1:Nt; % This is the case only for Tzer Han's data 
storeResult(numel(ct_fArr)) = struct();

% if ~exist('figSave', 'dir')
%    mkdir('figSave')
% end

% progressbar
   
f= figure('color','w','Units','normalized','OuterPosition',[0,0,1,1]);
for kk = 1:numel(ct_fArr)
    
    ct_f = ct_fArr(kk); ct_i = 1; 
    
    dt = 1*10^(-1); tau0_Advect = 0; tauf_Advect = time(ct_f)-time(ct_i); 
    xqBck = x{ct_f}; yqBck = y{ct_f}; zqBck = z{ct_f}; % Initial conditions for Backward FTLE 
    xqFor = x{ct_i}; yqFor = y{ct_i}; zqFor = z{ct_i}; %Initial condition for forward Advection- 2Check bFTLE results
    
    tf_Absolute = time(ct_f);
    tauArr_Advect = tau0_Advect:dt:tauf_Advect;
    tauSave_Advect = tauArr_Advect(1:round(numel(tauArr_Advect)/numel(tauArr_Advect)):end); %Can't save all the data, will save at selective times
    
    if isStatic == 1
        [xt_Advect,yt_Advect,zt_Advect] = AdvctBckStatic(cpu_num,tf_Absolute,tauSave_Advect,tauArr_Advect,xqBck,yqBck,zqBck,time,v,inputMesh);
    else
        [xt_Advect,yt_Advect,zt_Advect] = AdvctBck(cpu_num,tf_Absolute,tauSave_Advect,tauArr_Advect,xqBck,yqBck,zqBck,time,v,x,y,z,TrianT);
    end 
    
    % % Advect the grid forward
    
    t0_Advect_for = time(ct_i); tf_Advect_for = time(ct_f);
    tArr_Advect_for = t0_Advect_for:dt:tf_Advect_for;
    tSave_Advect_for = tArr_Advect_for(1:round(numel(tArr_Advect_for)/numel(tArr_Advect_for)):end); %Can't save all the data, will save at selective times
    
    if isStatic == 1
        [xt_Advect_for,yt_Advect_for,zt_Advect_for] = AdvctForStatic(cpu_num,tSave_Advect_for,tArr_Advect_for,xqFor,yqFor,zqFor,time,v,inputMesh);
    else
        [xt_Advect_for,yt_Advect_for,zt_Advect_for] = AdvctFor(cpu_num,tSave_Advect_for,tArr_Advect_for,xqFor,yqFor,zqFor,time,v,x,y,z,TrianT);
    end 
       
    %% Calculate the FTLE values 

    % Compute bw-FTLE
    x0 = squeeze(xt_Advect(1,:)); y0 = squeeze(yt_Advect(1,:)); z0 = squeeze(zt_Advect(1,:));
    xf = squeeze(xt_Advect(end,:)); yf = squeeze(yt_Advect(end,:)); zf = squeeze(zt_Advect(end,:));
    Tri_pts_Uni = TrianT{ct_f}; % The triangulation of the mesh at tf
    xfData = squeeze(x{1}); yfData = squeeze(y{1}); zfData = squeeze(z{1}); Trif_Data = squeeze(TrianT{1});
    
    [FTLE_vert,detdelF_vert,eigVec] = FTLEComputev3(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_pts_Uni,Trif_Data);

    storeResult(kk).x0Back = x0; storeResult(kk).y0Back = y0; storeResult(kk).z0Back = z0;
    storeResult(kk).xfBack = xf; storeResult(kk).yfBack = yf; storeResult(kk).zfBack = zf;
    storeResult(kk).xfDataBack = xfData; storeResult(kk).yfDataBack = yfData; storeResult(kk).zfDataBack = zfData;
    storeResult(kk).Tri_pts_UniBack = Tri_pts_Uni; storeResult(kk).Trif_DataBack = Trif_Data;
    storeResult(kk).FTLEBack = FTLE_vert; storeResult(kk).detdelFBack = detdelF_vert; storeResult(kk).eigVecBack = eigVec; 
    
    subplot(2,2,1)
    trisurf(Tri_pts_Uni,x0,y0,z0,log(detdelF_vert),'FaceAlpha',1,'Edgecolor','none');
    hold on; cmocean('balance','pivot',0); scatter3(squeeze(xt_Advect_for(end,:)),squeeze(yt_Advect_for(end,:)),squeeze(zt_Advect_for(end,:)),10,'r','filled')
    colorbar; axis equal; shading interp;  
    title('abs(detdelF)-bw')

    subplot(2,2,2)
    trisurf(Tri_pts_Uni,x0,y0,z0,log(FTLE_vert),'FaceAlpha',1,'Edgecolor','none'); hold on 
    quiver3(x0,y0,z0,eigVec(1,:),eigVec(2,:),eigVec(3,:),'r','ShowArrowHead','off'); hold off
    % hold on; scatter3(squeeze(xt_Advect_for(end,:)),squeeze(yt_Advect_for(end,:)),squeeze(zt_Advect_for(end,:)),'r','filled')
    colorbar; axis equal; shading interp;  
    title('bwFTLE')
    
    % Compute FW-FTLE
    x0 = squeeze(xt_Advect_for(1,:)); y0 = squeeze(yt_Advect_for(1,:)); z0 = squeeze(zt_Advect_for(1,:));
    xf = squeeze(xt_Advect_for(end,:)); yf = squeeze(yt_Advect_for(end,:)); zf = squeeze(zt_Advect_for(end,:));
    Tri_pts_Uni = TrianT{ct_i}; xfData = squeeze(x{ct_f}); yfData = squeeze(y{ct_f}); zfData = squeeze(z{ct_f}); Trif_Data = squeeze(TrianT{ct_f});
    
    [FTLE_vert,detdelF_vert,eigVec] = FTLEComputev3(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_pts_Uni,Trif_Data);
    
    storeResult(kk).x0For = x0; storeResult(kk).y0For = y0; storeResult(kk).z0For = z0;
    storeResult(kk).xfFor = xf; storeResult(kk).yfFor = yf; storeResult(kk).zfFor = zf;
    storeResult(kk).xfDataFor = xfData; storeResult(kk).yfDataFor = yfData; storeResult(kk).zfDataFor = zfData;
    storeResult(kk).Tri_pts_UniFor = Tri_pts_Uni; storeResult(kk).Trif_DataFor = Trif_Data;
    storeResult(kk).FTLEFor = FTLE_vert; storeResult(kk).detdelFFor = detdelF_vert; storeResult(kk).eigVecFor = eigVec;
    storeResult(kk).t0 = time(ct_i); storeResult(kk).tf = time(ct_f);

    subplot(2,2,3)
    trisurf(Tri_pts_Uni,x0,y0,z0,log(detdelF_vert),'FaceAlpha',1,'Edgecolor','none');
    cmocean('balance','pivot',0);colorbar; axis equal; shading interp;  
    title('abs(detdelF)-fw')

    subplot(2,2,4)
    trisurf(Tri_pts_Uni,x0,y0,z0,log(FTLE_vert),'FaceAlpha',1,'Edgecolor','none');hold on 
    quiver3(x0,y0,z0,eigVec(1,:),eigVec(2,:),eigVec(3,:),'r','ShowArrowHead','off'); hold off
    colorbar; axis equal; shading interp;  
    title('fwFTLE')

    % sgtitle(sprintf('Time = %d',ct_f))
    disp(ct_f)
    % progressbar(kk/numel(ct_fArr))
end 
save('./structData','storeResult')
