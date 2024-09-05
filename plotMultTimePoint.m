% Plot the result of the spheroids for multiple time points 
clear; clc; close all; 
load('.\structData.mat')
Nt = size(storeResult,2);

if ~exist('figSave', 'dir')
   mkdir('figSave')
end

camView = [10,-71];camvaAmp = 9;
numNNSmooth = 10;
%% Plot the results of the multTime FTLE calc - backward and forward FTLE 

if ~exist('test', 'dir')
    mkdir('./test')
else
    rmdir('./test','s')
    mkdir('./test')
end 

% Create mask for the deformation Grid at init time 
load('./formatData/dataHeart.mat')

createMaskAtInitTime = [storeResult(1).x0For',storeResult(1).y0For',storeResult(1).z0For'];
connecAtInitTime = storeResult(1).Tri_pts_UniFor;
x0DefCenter = 436; y0DefCenter = 487; z0DefCenter = 107;
IdxDefGridVertices = knnsearch(createMaskAtInitTime,[x0DefCenter,y0DefCenter,z0DefCenter] ...
    ,'K',1000);
% trisurf(connecAtInitTime,createMaskAtInitTime(:,1),createMaskAtInitTime(:,2),createMaskAtInitTime(:,3))
% hold on; scatter3(createMaskAtInitTime(IdxDefGridVertices,1),createMaskAtInitTime(IdxDefGridVertices,2), ...
%     createMaskAtInitTime(IdxDefGridVertices,3),'r','filled')

clc; close all; VideoFileName=strcat(cd,sprintf('/multTimeEverything'));  % How much to coarsen velocity field
writerObj = VideoWriter(VideoFileName,'MPEG-4');writerObj.FrameRate = 5;  writerObj.Quality = 100;
writerObj.FileFormat; open(writerObj); count = 0; fntSz = 30; coarse = 1; 

f = figure('color','w','Units','normalized','OuterPosition',[0.0100 0.0178 0.9799 0.4744]);

for kk = 1:Nt
    v1Data = squeeze(v{1,kk+1}); v2Data = squeeze(v{2,kk+1}); v3Data = squeeze(v{3,kk+1});
    xData = squeeze(x{kk+1}); yData = squeeze(y{kk+1}); zData = squeeze(z{kk+1});
    TriData = TrianT{kk+1}; vMag_Data = sqrt(v1Data.^2+v2Data.^2+v3Data.^2);
    
    subplot(1,3,1)
    p = trisurf(TriData,xData,yData,zData,vMag_Data,'FaceAlpha',1,'Edgecolor','none'); hold on 
    p.SpecularStrength = 0.2; p.AmbientStrength = 0.3; c = colorbar; shading interp
    quiver3(xData(1:coarse:end),yData(1:coarse:end),zData(1:coarse:end),...
        v1Data(1:coarse:end),...
        v2Data(1:coarse:end),...
        v3Data(1:coarse:end),...
        'LineWidth',1,'Color','r','MaxHeadSize',20)
    c.FontSize = fntSz;view(camView); camva(camvaAmp); axis equal off; camlight; clim([0,26])
    set(gcf,'color','w'); title(sprintf('$$ \\mathbf{v}(\\mathbf{x},t = %d) $$',kk),'Interpreter','latex','FontSize',fntSz);
    hold off

    %% Plot the FTLE 
    count = count+1;
    x0 = storeResult(kk).x0Back; y0 = storeResult(kk).y0Back; z0 = storeResult(kk).z0Back;
    Tri_pts_Uni = storeResult(kk).Tri_pts_UniBack;
    FTLE_vert = storeResult(kk).FTLEBack; detdelF_vert = storeResult(kk).detdelFBack; eigVec = storeResult(kk).eigVecBack;
    ct_f = storeResult(kk).tf;ct_0 = storeResult(kk).t0;   
    logicArr = FTLE_vert>7;
    
    [Idx,~] = knnsearch([x0',y0',z0'],[x0',y0',z0'],"K",numNNSmooth);
    FTLE_vertSmooth = mean(FTLE_vert(Idx),2);
    FTLE_vert(logicArr) = FTLE_vertSmooth(logicArr); % Replace the high value regions 

    % [plotx,ploty,plotz] = findOuterBoundaryHeart(x0,y0,z0,Tri_pts_Uni);
    % ^Ready- But shows the backbone too (How to fix>?) 
    
    subplot(1,3,2)
    p = trisurf(Tri_pts_Uni,x0,y0,z0,FTLE_vert,'FaceAlpha',1,'Edgecolor','black'); hold on 
    p.SpecularStrength = 0.2; p.AmbientStrength = 0.3;
    % plot3(plotx,ploty,plotz,'k','LineWidth',2); hold on; % Plots the ends and backbone
    % quiver3(x0(1:crField:end),y0(1:crField:end),z0(1:crField:end),eigVec(1,1:crField:end),eigVec(2,1:crField:end),eigVec(3,1:crField:end),'b','ShowArrowHead','off','LineWidth',1); hold on
    c = colorbar; axis equal off;  shading interp; camva(camvaAmp);
    c.FontSize = fntSz; view(camView); camlight; clim([0,7])
    % title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}(\\mathbf{x}_t), \\mathbf{\\xi}_{%.1f}^{%.1f}$$',ct_f-ct_0,0,ct_f-ct_0,0),'Interpreter','latex','FontSize',fntSz)
    title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}(\\mathbf{x}_t)$$',ct_f-ct_0,0),'Interpreter','latex','FontSize',fntSz)
    hold off
    % %% Plot the deformation Grid
    % reducFac = 0.5;
    % lineSentivity = 0.05;
    % 
    % % % Lot of projections to do - Prepare the Mesh Object
    % inputMesh = struct();
    % xMesh = squeeze(storeResult(kk).x0Back)'; yMesh = squeeze(storeResult(kk).y0Back)'; zMesh = squeeze(storeResult(kk).z0Back)';
    % inputMesh.faces = storeResult(kk).Tri_pts_UniBack; inputMesh.nodes = [squeeze(xMesh), squeeze(yMesh),squeeze(zMesh)];
    % 
    % resultDefGrid = defGridCalcHeart(IdxDefGridVertices,inputMesh, ...
    %     storeResult(kk).x0For, ...
    %     storeResult(kk).y0For,storeResult(kk).z0For, ...
    %     storeResult(kk).xfFor,storeResult(kk).yfFor,storeResult(kk).zfFor, ...
    %     storeResult(kk).Tri_pts_UniFor,reducFac,lineSentivity);
    % 
    % xdefGrid = squeeze(resultDefGrid(:,1,:))';
    % ydefGrid = squeeze(resultDefGrid(:,2,:))';
    % zdefGrid = squeeze(resultDefGrid(:,3,:))';
    % 
    % xfGrid = storeResult(kk).xfFor;
    % yfGrid = storeResult(kk).yfFor;
    % zfGrid = storeResult(kk).zfFor;
    % 
    % subplot(1,3,2)
    % p = trisurf(Tri_pts_Uni,x0,y0,z0,FTLE_vert,'FaceAlpha',1,'Edgecolor','black'); hold on 
    % p.SpecularStrength = 0.2; p.AmbientStrength = 0.3;
    % % scatter3(xfGrid(IdxDefGridVertices),yfGrid(IdxDefGridVertices),zfGrid(IdxDefGridVertices),'r','filled')
    % % plot3(plotx,ploty,plotz,'k','LineWidth',2); hold on; % Plots the ends and backbone
    % plot3(xdefGrid,ydefGrid,zdefGrid,'r','LineWidth',1); hold on 
    % % quiver3(x0(1:crField:end),y0(1:crField:end),z0(1:crField:end),eigVec(1,1:crField:end),eigVec(2,1:crField:end),eigVec(3,1:crField:end),'b','ShowArrowHead','off','LineWidth',1); hold on
    % c = colorbar; axis equal off;  shading interp; camva(camvaAmp);
    % c.FontSize = fntSz; view(camView); camlight;
    % % title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}(\\mathbf{x}_t), \\mathbf{\\xi}_{%.1f}^{%.1f}$$',ct_f-ct_0,0,ct_f-ct_0,0),'Interpreter','latex','FontSize',fntSz)
    % title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}(\\mathbf{x}_t)$$',ct_f-ct_0,0),'Interpreter','latex','FontSize',fntSz)
    % hold off
    
    %% Forward FTLE 
    x0 = storeResult(kk).x0For; y0 =  storeResult(kk).y0For;z0 =  storeResult(kk).z0For;
    xf = storeResult(kk).xfFor; yf = storeResult(kk).yfFor; zf = storeResult(kk).zfFor;
    xfData = storeResult(kk).xfDataFor;yfData =  storeResult(kk).yfDataFor;zfData =  storeResult(kk).zfDataFor;
    Tri_pts_Uni = storeResult(kk).Tri_pts_UniFor;Trif_Data =  storeResult(kk).Trif_DataFor;
    FTLE_vert = storeResult(kk).FTLEFor; detdelF_vert = storeResult(kk).detdelFFor;eigVec = storeResult(kk).eigVecFor;
    
    % [plotx,ploty,plotz] = findOuterBoundaryHeart(x0,y0,z0,Tri_pts_Uni);

    f1 = subplot(1,3,3);
    p = trisurf(Tri_pts_Uni,x0,y0,z0,log(FTLE_vert),'FaceAlpha',1,'Edgecolor','none');hold on 
    p.SpecularStrength = 0.2; p.AmbientStrength = 0.3;
    % plot3(plotx,ploty,plotz,'k','LineWidth',2); hold on; % Plots the ends and backbone
    % quiver3(x0(1:crField:end),y0(1:crField:end),z0(1:crField:end),eigVec(1,1:crField:end),eigVec(2,1:crField:end),eigVec(3,1:crField:end),'b','ShowArrowHead','off','LineWidth',1); hold off
    c = colorbar; axis equal off; shading interp; camva(camvaAmp);
    c.FontSize = fntSz; view(camView); camlight; clim([-2.5,3])
    % title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}(\\mathbf{x}_t), \\mathbf{\\xi}_{%.1f}^{%.1f}$$',0,ct_f-ct_0,0,ct_f-ct_0),'Interpreter','latex','FontSize',fntSz)
    title(sprintf('$$ \\Lambda_{%.1f}^{%.1f}(\\mathbf{x}_t)$$',0,ct_f-ct_0),'Interpreter','latex','FontSize',fntSz)
    hold off

    fname = ['.\test\test' num2str(count)]; % full name of image
    print('-djpeg','-r300',fname)     % save image with '-r200' resolution
    I = imread([fname '.jpg']);       % read saved image
    frame = im2frame(I);              % convert image to frame
    writeVideo(writerObj,frame);

end 
close(writerObj);
