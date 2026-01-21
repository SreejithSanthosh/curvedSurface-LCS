function [face_mean_nodes,face_normals]=getFaceCenterAndNormals_new(faces,nodes)
% determines the incenter of a face and the normal direction to the face    
TR = triangulation(faces,nodes);
face_mean_nodes = incenter(TR);
face_normals = faceNormal(TR);
end