function [v1q,v2q,v3q] = interpOnMesh(Triangul,xData,yData,zData,v1Data,v2Data,v3Data,xq,yq,zq)
% Intialize the variables to store 
v1q = 0*xq; v2q = 0*xq; v3q = 0*xq;

%Making the triangulation for the mesh
FV.nodes = [squeeze(xData), squeeze(yData),squeeze(zData)]; %Faces
FV.faces = Triangul; points = [xq,yq,zq];
[~,baryCoord,faceId] = fastPoint2TriMeshSRJv1(FV,points,0);

for i = 1:size(baryCoord,1)
    faceId_i = faceId(i); IdxFaceVert = Triangul(faceId_i,:);
    bar1 = baryCoord(i,1); bar2 = baryCoord(i,2);

    v1q(i) = (1-bar1-bar2)*v1Data(IdxFaceVert(1))+bar1*v1Data(IdxFaceVert(2))+bar2*v1Data(IdxFaceVert(3));
    v2q(i) = (1-bar1-bar2)*v2Data(IdxFaceVert(1))+bar1*v2Data(IdxFaceVert(2))+bar2*v2Data(IdxFaceVert(3));
    v3q(i) = (1-bar1-bar2)*v3Data(IdxFaceVert(1))+bar1*v3Data(IdxFaceVert(2))+bar2*v3Data(IdxFaceVert(3));
end 


end 