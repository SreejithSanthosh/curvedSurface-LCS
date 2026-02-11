function [s2,s1,e2,e1] = compute_strain(mesh_r0,mesh_F0,mesh_v0,delta)
% This code also computes the eigenvectors and eigenvalus of  the
% strain-rate tensor to identify the OECS

% mesh_r0 : (Nv_0 x 3) double array
% mesh_F0 : (Nf_0 x 3) double array
% mesh_v0 : (Nv_0 x 3) double array
% delta : scalar double : spatial smoothing

G = meshToGraph(mesh_r0, mesh_F0); % Graph for finding nearby tracers

Nq = size(mesh_r0,1); % Number of points being advected
s2 = nan(Nq, 1); % Largest eigenvalue of B'*B
s1 = nan(Nq, 1); % Smallest eigenvalue of B'*B
e1 = nan(Nq, 3); % Initialize eigenvector array at x0
e2 = nan(Nq, 3); % Initialize eigenvector array at xf

TR_0 = triangulation(mesh_F0,mesh_r0);
VN_0q = vertexNormal(TR_0); % Compute normal

for i = 1:Nq
     % Computes the nearest nodes within this distance
    neighb = nearest(G,i,delta);
    
    % Consitruct the matrix 
    X = mesh_r0(neighb,:)-mesh_r0(i,:); Y = mesh_v0(neighb,:)-mesh_v0(i,:);
    
    % Remove the normal component 
    n0 = repmat(VN_0q(i,:),size(X,1),1); 

    n0comp = sum(X.*n0,2);
    X = X - repmat(n0comp,1,3).*n0;
    Y = Y - repmat(n0comp,1,3).*n0;
    
    warning('off','MATLAB:rankDeficientMatrix')
    grad_v = Y'/X'; % (3 x 3) double  
    warning('on','MATLAB:rankDeficientMatrix')
    
    % Create tangent vector basis at t0
    n0q = VN_0q(i,:);
    if abs(dot(n0q,[1,0,0]))>0.9
        arb_vec = [0,1,0];
    else
        arb_vec = [1,0,0];
    end
    zeta1_0 = cross(n0q,arb_vec);
    zeta2_0 = cross(n0q,zeta1_0);
    zeta1_0 = zeta1_0./norm(zeta1_0);
    zeta2_0 = zeta2_0./norm(zeta2_0);

    E = nan(2,2); % E is gradV projected onto the tangent vector basis 
    E(1,1) = zeta1_0*grad_v*zeta1_0'; E(1,2) = zeta1_0*grad_v*zeta2_0';
    E(2,1) = zeta2_0*grad_v*zeta1_0'; E(2,2) = zeta2_0*grad_v*zeta2_0';
    
    D = 0.5*(E+E');
    [eV0,eigMatr] = eig(D);
    [eiglist,I] = sort(diag(eigMatr)); 

    eV0 = eV0(:,I); 
    V_large0 = eV0(1,2)*zeta1_0+eV0(2,2)*zeta2_0;
    V_small0 = eV0(1,1)*zeta1_0+eV0(2,1)*zeta2_0;
       
    % Compute the L2 field
    s2(i) = max(eiglist); % largest eig. value
    s1(i) = min(eiglist); % smallest eig. value
    e2(i,:) = V_large0; % Store the eigVec of largest eig val. at x0
    e1(i,:) = V_small0; % Store the eigVec of smallest eig val. at x0

end 



