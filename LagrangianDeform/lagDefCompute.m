function [FTLE_vert,detdelF_vert,eig_vert] = lagDefCompute(x0,y0,z0,xf,yf,zf,xfData,yfData,zfData,Tri_q0,Tri_qf,tf,ti)

% (x0,y0,z0) is the intial tracer location 
% (xf,yf,zf) is the final tracer location 
% (xfData, yfData, zfData) is the given mesh at the final location 
% (Tri_q0,Tri_qf) is the mesh at t = t0 and t = tf

Nvert = size(x0,2); 
FTLE_vert = nan([Nvert,1]); detdelF_vert = nan([Nvert,1]);
eig_vert = nan([3,Nvert]); vertices = [x0',y0',z0'];
mesh = surfaceMesh(vertices,Tri_q0); computeNormals(mesh,'vertex');

% Find normal at final time 
nvecfArr = normalAtPointGivenMesh(xf',yf',zf',xfData,yfData,zfData,Tri_qf);

for i = 1:Nvert
    x0vert = x0(i); y0vert = y0(i); z0vert = z0(i);
    xfvert = xf(i); yfvert = yf(i); zfvert = zf(i);
    
    % Find list of all neighbors     
    faceWithVertex = logical(sum(Tri_q0 == i,2)); idXUnique = unique(Tri_q0(faceWithVertex,:));
    idXUnique = idXUnique(idXUnique ~= i); zetaArr = [x0(idXUnique)-x0vert;y0(idXUnique)-y0vert;z0(idXUnique)-z0vert];
    betaArr = [xf(idXUnique)-xfvert;yf(idXUnique)-yfvert;zf(idXUnique)-zfvert];
    
    % % Compute the normal vector at intial time 
    nVec0 = mesh.VertexNormals(i,:)';

    % Find normal at final time 
    % nvecf = normalAtPointGivenMesh(xfvert,yfvert,zfvert,xfData,yfData,zfData,Tri_qf);
    nvecf = nvecfArr(:,i);

    % Construct the vectors to get rid of the normal direction 

    for jj = 1: size(zetaArr,2) % Remove normal from the intial vector configuration 
        zeta_jj = zetaArr(:,jj);
        zeta_jj = zeta_jj-(nVec0'*zeta_jj).*nVec0;

        zetaArr(:,jj) = zeta_jj;
    end

    for jj = 1: size(betaArr,2) % Remove normal from the final vector configuration 
        beta_jj = betaArr(:,jj);
        beta_jj = beta_jj-(nvecf'*beta_jj).*nvecf;

        betaArr(:,jj) = beta_jj;
    end
    
    % Rotate the vectors to normal along z
    % Initial configuration
    if ~all(nVec0 ==[0;0;1])
        rotAxis = cross(nVec0,[0;0;1]);
    else 
        rotAxis = [0;0;1];
    end 

    rotAng0 = atan(sqrt(nVec0(1)^2+nVec0(2)^2)/nVec0(3));
    if abs(rotAng0) > exp(-6) 
        M = AxelRot(rad2deg(rotAng0),rotAxis',[0,0,0]); rotMatr0 = M(1:3,1:3); 
    else
        rotMatr0 = eye(3);
    end
      
    % Final configuration
    if ~all( nvecf==[0;0;1] )
        rotAxis = cross(nvecf,[0;0;1]);
    else 
        rotAxis = [0;0;1];
    end

    rotAngf = atan(sqrt(nvecf(1)^2+nvecf(2)^2)/nvecf(3));   
    if abs(rotAngf) > exp(-6) 
        M = AxelRot(rad2deg(rotAngf),rotAxis',[0,0,0]); rotMatrf = M(1:3,1:3);
    else
        rotMatrf = eye(3);
    end 
    
    zetaArr = rotMatr0*zetaArr; betaArr = rotMatrf*betaArr;
    
    X = zetaArr(1:2,:); Y = betaArr(1:2,:);
    regulFac = 10^(-7); delF = (Y*X'+regulFac*size(X,2)*eye(2))*inv(X*X'+regulFac*size(X,2)*eye(2));
    cauchMatr = delF'*delF; [V,e] = eig(cauchMatr);
    FTLE_vert(i) = log(sqrt(max(diag(e))))/abs(tf-ti); 
    detdelF_vert(i) = log(abs(det(delF)))/abs(tf-ti);

    [~,I] = sort(diag(e)); V = V(:,I); eigMax = zeros(3,1); eigMax(1:2) = V(:,2);
    eig_vert(:,i) = rotMatr0'*eigMax;

end 


end 