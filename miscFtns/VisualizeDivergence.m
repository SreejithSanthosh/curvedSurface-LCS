% We visaulize the divergence of the velocity fields 
clear; clc; close all;
load('../Data/processedS1.mat')
load('../Data/S1_OmegaMat_v2.mat')
Nt = size(TrianT,1);

%% Make figure for the paper 
clc; close all;

fntSz = 30; camView = [60,7.2];
f = figure('color','w','Units','normalized','OuterPosition',[0.0031 0.0056 0.9938 0.5500]);

i = 40; % Select the time point at which you want to visualize
omega_i = OmegaMat(:,i); omega_i = omega_i./norm(omega_i);
rOmega = 1.2; xOmega = -rOmega.*[-omega_i(1),omega_i(1)];
yOmega = -rOmega.*[-omega_i(2),omega_i(2)]; zOmega = -rOmega.*[-omega_i(3),omega_i(3)];

vi_x = v{1,i}; vi_y = v{2,i}; vi_z = v{3,i}; 
xi = x{i}; yi = y{i}; zi = z{i}; triI = TrianT{i};
Nfaces = size(triI,1); div_i = zeros([Nfaces,1]);
TR = triangulation(triI,xi,yi,zi); C = incenter(TR);
xinCen = C(:,1); yinCen = C(:,2); zinCen = C(:,3);
convCen = convhulln([xinCen,yinCen,zinCen]);

for j = 1:Nfaces
    idx1 = triI(j,1); idx2 = triI(j,2); idx3 = triI(j,3);

    zeta1 = [xi(idx2)-xi(idx1),yi(idx2)-yi(idx1),zi(idx2)-zi(idx1)];
    zeta2 = [xi(idx3)-xi(idx1),yi(idx3)-yi(idx1),zi(idx3)-zi(idx1)];
    
    del_zeta1 = [vi_x(idx2)-vi_x(idx1),vi_y(idx2)-vi_y(idx1),vi_z(idx2)-vi_z(idx1)];
    del_zeta2 = [vi_x(idx3)-vi_x(idx1),vi_y(idx3)-vi_y(idx1),vi_z(idx3)-vi_z(idx1)];
    
    a = cross(zeta1,zeta2); b = cross(del_zeta1,zeta2)+cross(zeta1,del_zeta2);
    div_i(j) = dot(a,b)./norm(a)^2;       
end 

subplot(1,2,1)
trisurf(convCen,xinCen,yinCen,zinCen,div_i,'Edgecolor','none'); axis equal; xlabel('x'); ylabel('y'); zlabel('z');
cmocean('balance','pivot',0); 
% title(sprintf('time = %d, meanDiv = %f',i, mean(div_i))); 
c = colorbar; axis equal off; shading interp; camva(9); c.FontSize = fntSz; view(camView);
hold on; plot3(xOmega,yOmega,zOmega,'k','LineWidth',2); hold off;
title(sprintf('$$ \\mathbf{\\nabla}\\cdot\\mathbf{v}(\\mathbf{x},t = %d) $$',i),'Interpreter','latex','FontSize',fntSz)

subplot(1,2,2)
trisurf(convCen,xinCen,yinCen,zinCen,div_i,'Edgecolor','none'); xlabel('x'); ylabel('y'); zlabel('z'); axis equal
cmocean('balance','pivot',0); 
% title(sprintf('time = %d',i)); 
c = colorbar; axis equal off; shading interp; camva(9); c.FontSize = fntSz; view(camView);
hold on; plot3(xOmega,yOmega,zOmega,'k','LineWidth',2); hold off;
title(sprintf('$$ <\\mathbf{\\nabla}\\cdot\\mathbf{v}(\\mathbf{x},t = %d)> = %.4f$$',i,mean(div_i)),'Interpreter','latex','FontSize',fntSz)

export_fig(sprintf('./divFigPaper%d.png',i),'-m2');
%% Plot the divergence 
clear mov; clear writerObj; clc; close all;
VideoFileName=strcat('./divg'); writerObj = VideoWriter(VideoFileName,'MPEG-4');
writerObj.FrameRate = 3;  writerObj.Quality = 100; count = 0; open(writerObj); fntSz = 30; camView = [60,7.2];

f = figure('color','w','Units','normalized','OuterPosition',[0.0031 0.0056 0.9938 0.5500]);
for i = 1:Nt
    count = count+1; omega_i = OmegaMat(:,i); omega_i = omega_i./norm(omega_i);
    rOmega = 1.2; xOmega = -rOmega.*[-omega_i(1),omega_i(1)];
    yOmega = -rOmega.*[-omega_i(2),omega_i(2)]; zOmega = -rOmega.*[-omega_i(3),omega_i(3)];

    vi_x = v{1,i}; vi_y = v{2,i}; vi_z = v{3,i}; 
    xi = x{i}; yi = y{i}; zi = z{i}; triI = TrianT{i};
    Nfaces = size(triI,1); div_i = zeros([Nfaces,1]);
    TR = triangulation(triI,xi,yi,zi); C = incenter(TR);
    xinCen = C(:,1); yinCen = C(:,2); zinCen = C(:,3);
    convCen = convhulln([xinCen,yinCen,zinCen]);

    for j = 1:Nfaces
        idx1 = triI(j,1); idx2 = triI(j,2); idx3 = triI(j,3);

        zeta1 = [xi(idx2)-xi(idx1),yi(idx2)-yi(idx1),zi(idx2)-zi(idx1)];
        zeta2 = [xi(idx3)-xi(idx1),yi(idx3)-yi(idx1),zi(idx3)-zi(idx1)];
        
        del_zeta1 = [vi_x(idx2)-vi_x(idx1),vi_y(idx2)-vi_y(idx1),vi_z(idx2)-vi_z(idx1)];
        del_zeta2 = [vi_x(idx3)-vi_x(idx1),vi_y(idx3)-vi_y(idx1),vi_z(idx3)-vi_z(idx1)];
        
        a = cross(zeta1,zeta2); b = cross(del_zeta1,zeta2)+cross(zeta1,del_zeta2);
        div_i(j) = dot(a,b)./norm(a)^2;       
    end 
    
    subplot(1,2,1)
    trisurf(convCen,xinCen,yinCen,zinCen,div_i,'Edgecolor','none'); axis equal; xlabel('x'); ylabel('y'); zlabel('z');
    cmocean('balance','pivot',0); 
    % title(sprintf('time = %d, meanDiv = %f',i, mean(div_i))); 
    c = colorbar; axis equal off; shading interp; camva(9); c.FontSize = fntSz; view(camView);
    hold on; plot3(xOmega,yOmega,zOmega,'k','LineWidth',2); hold off;
    title(sprintf('$$ \\mathbf{\\nabla}\\cdot\\mathbf{v}(\\mathbf{x},t = %d) $$',i),'Interpreter','latex','FontSize',fntSz)
    
    subplot(1,2,2)
    trisurf(convCen,xinCen,yinCen,zinCen,div_i,'Edgecolor','none'); xlabel('x'); ylabel('y'); zlabel('z'); axis equal
    cmocean('balance','pivot',0); 
    % title(sprintf('time = %d',i)); 
    c = colorbar; axis equal off; shading interp; camva(9); c.FontSize = fntSz; view(camView);
    hold on; plot3(xOmega,yOmega,zOmega,'k','LineWidth',2); hold off;
    title(sprintf('$$ <\\mathbf{\\nabla}\\cdot\\mathbf{v}(\\mathbf{x},t = %d)> = %.4f$$',i,mean(div_i)),'Interpreter','latex','FontSize',fntSz)

    mov(count) = getframe(gcf);


end 

writeVideo(writerObj,mov); close(writerObj); 