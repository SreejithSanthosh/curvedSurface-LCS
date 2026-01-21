function vq = interp_on_mesh(mesh_F,mesh_r,mesh_v,rq)
arguments
    mesh_F (:,3) double
    mesh_r (:,3) double
    mesh_v  
    rq (:,3) double
end 

% The mesh_v = (Nv x ** ): can be any dimension 
%Making the triangulation for the mesh
FV.nodes = mesh_r; % (Nf x 3)
FV.faces = mesh_F; % (Nv x 3)
points = rq; % (Nq x 3)

[~,baryCoord,faceId] = fastPoint2TriMeshModif(FV,points,0);

% 1. Extract the indices for all vertices of all relevant faces at once
% mesh_F(faceId, :) returns an Nq x 3 matrix of vertex indices
IdxFaceVerts = mesh_F(faceId, :);

% 2. Get the specific vertex indices (q1, q2, q3) for every query point
q1 = IdxFaceVerts(:, 1);
q2 = IdxFaceVerts(:, 2);
q3 = IdxFaceVerts(:, 3);

% 3. Calculate the weights for each vertex
% bar1 and bar2 are Nq x 1 columns; w1 becomes Nq x 1
w1 = 1 - baryCoord(:, 1) - baryCoord(:, 2);
w2 = baryCoord(:, 1);
w3 = baryCoord(:, 2);

% 4. Compute the interpolated positions (vq)
% We use the weights to scale the rows of mesh_v indexed by q1, q2, and q3
vq = w1 .* mesh_v(q1, :) + ...
     w2 .* mesh_v(q2, :) + ...
     w3 .* mesh_v(q3, :);

% for i = 1:Nq
%     faceId_i = faceId(i); IdxFaceVert = mesh_F(faceId_i,:);
%     bar1 = baryCoord(i,1); bar2 = baryCoord(i,2);
% 
% 
%      q1 = IdxFaceVert(1);
%      q2 = IdxFaceVert(2);
%      q3 = IdxFaceVert(3);
% 
%      vq(i,:) = (1-bar1-bar2)*mesh_v(q1,:)+bar1*mesh_v(q2,:)+bar2*mesh_v(q3,:);
% end 

end 