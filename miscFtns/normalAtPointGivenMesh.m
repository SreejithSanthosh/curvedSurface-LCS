function nvecf = normalAtPointGivenMesh(xq,yq,zq,xData,yData,zData,Triangul)
% nvecf : (3, #points) 
nvecf = zeros([3,numel(xq)]);

%Making the triangulation for the mesh
FV.nodes = [squeeze(xData), squeeze(yData),squeeze(zData)]; %Faces
FV.faces = Triangul; points = [xq,yq,zq];

[~,~,faceId] = fastPoint2TriMeshSRJv1(FV,points,0);

for i = 1:numel(xq)
    faceId_i = faceId(i); idXMesh = Triangul(faceId_i,:);

    coord1 = [xData(idXMesh(1)),yData(idXMesh(1)),zData(idXMesh(1))];
    coord2 = [xData(idXMesh(2)),yData(idXMesh(2)),zData(idXMesh(2))];
    coord3 = [xData(idXMesh(3)),yData(idXMesh(3)),zData(idXMesh(3))];

    zeta1 = coord2-coord1; zeta2 = coord3-coord1;
    nvecf(:,i) = cross(zeta1,zeta2)./norm(cross(zeta1,zeta2));
end 

end 