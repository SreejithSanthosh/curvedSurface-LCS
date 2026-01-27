function [L2,L1,V0,Vf] = compute_deform(mesh_r0,mesh_rf,mesh_F0,mesh_Ff,rf,delta)
% This code also computes the eigenvectors
% mesh_r0 : (Nv_0 x 3) double array
% mesh_rf : (Nv_f x 3) double array
% mesh_F0 : (Nf_0 x 3) double array
% mesh_Ff : (Nf_f x 3) double array
% rf: (Nq x 3) double array
% delta : scalar double : spatial smoothing

r0 = mesh_r0; % Mesh at t0 provides init conditions for tracers
G = meshToGraph(r0, mesh_F0); % Graph for finding nearby tracers

% Project final position of tracer particle 
% to the final time position-mesh (For consistency)
mesh_struct_f.faces = mesh_Ff;
mesh_struct_f.nodes = mesh_rf;
rf = fastPoint2TriMeshModif(mesh_struct_f,rf);

Nq = size(r0,1); % Number of points being advected
L2 = nan(Nq, 1); % Initialize L2 array
L1 = nan(Nq, 1); % Initialize lowest eigenvalue 
V0 = nan(Nq, 3); % Initialize eigenvector array at x0
Vf = nan(Nq, 3); % Initialize eigenvector array at xf


TR_0 = triangulation(mesh_F0,mesh_r0);
TR_f = triangulation(mesh_Ff,mesh_rf);

VN_0 = vertexNormal(TR_0); % Compute normal
VN_f = vertexNormal(TR_f); % Compute normal
VN_0q = interp_on_mesh(mesh_F0,mesh_r0,VN_0,r0); % Interpolate t0
VN_fq = interp_on_mesh(mesh_Ff,mesh_rf,VN_f,rf); % Interpolate tf

for i = 1:Nq
     % Computes the nearest nodes within this distance
    neighb = nearest(G,i,delta);
    
    % Consitruct the matrix 
    X = r0(neighb,:)-r0(i,:); Y = rf(neighb,:)-rf(i,:);
    
    % Remove the normal component 
    n0 = repmat(VN_0q(i,:),size(X,1),1); 
    nf = repmat(VN_fq(i,:),size(X,1),1);
    
    n0comp = sum(X.*n0,2); nfcomp = sum(Y.*nf,2);
    X = X - repmat(n0comp,1,3).*n0;
    Y = Y - repmat(nfcomp,1,3).*nf;
    
    warning('off','MATLAB:rankDeficientMatrix')
    DF = Y'/X'; % (3 x 3) double  
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

    % Create tangent vector basis at tf
    nfq = VN_fq(i,:);
    if abs(dot(nfq,[1,0,0]))>0.9
        arb_vec = [0,1,0];
    else
        arb_vec = [1,0,0];
    end
    zeta1_f = cross(nfq,arb_vec);
    zeta2_f = cross(nfq,zeta1_f);
    zeta1_f = zeta1_f./norm(zeta1_f);
    zeta2_f = zeta2_f./norm(zeta2_f);

    B = nan(2,2); % B is DF projected onto the tangent vector basis 
    B(1,1) = zeta1_f*DF*zeta1_0'; B(1,2) = zeta1_f*DF*zeta2_0';
    B(2,1) = zeta2_f*DF*zeta1_0'; B(2,2) = zeta2_f*DF*zeta2_0';
    

    [eV0,D] = eig(B'*B);
    [eiglist,I] = sort(diag(D)); 
    eV0 = eV0(:,I); V_large0 = eV0(1,2)*zeta1_0+eV0(2,2)*zeta2_0;
    eVf = B*eV0; V_largef = eVf(1,2)*zeta1_f+eVf(2,2)*zeta2_f;
    
    % Compute the L2 field
    L2(i) = log(max(eiglist)); % L2 calculation
    L1(i) = log(min(eiglist)); % lowest eigenval. calculation
    V0(i,:) = V_large0; % Store the eigVec of largest eig val. at x0
    Vf(i,:) = V_largef; % Store the eigVec of largest eig val. at xf

end 


end 
