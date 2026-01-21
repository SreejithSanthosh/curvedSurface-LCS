function [FTLE,J] = compute_deformation_basic(mesh_r0,mesh_rf,mesh_F0,mesh_Ff,rf,delta)
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
FTLE = nan(Nq, 1); % Initialize FTLE array
J = nan(Nq, 1); % Initialize Area change array

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
    

    [~,D] = eig(B'*B);
    eiglist = sqrt(sort(diag(D)));

    % Compute the FTLE field
    J(i) = log(abs(eiglist(1)*eiglist(2))); % Jacobian determinant
    FTLE(i) = log(max(eiglist)); % FTLE calculation
end 


end 
