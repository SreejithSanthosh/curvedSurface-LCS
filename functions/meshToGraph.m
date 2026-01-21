function G = meshToGraph(vertices, faces)
    % Vertices: N x 3 matrix of coordinates
    % Faces:    M x 3 (triangular) or M x 4 (quad) matrix of vertex indices

    % 1. Extract all edges from the faces
    % For a triangular mesh, edges are (1,2), (2,3), and (3,1)
    if size(faces, 2) == 3
        edgePairs = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
    elseif size(faces, 2) == 4
        % For quad meshes
        edgePairs = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 4]); faces(:, [4, 1])];
    else
        error('Unsupported face format. Faces should be Nx3 or Nx4.');
    end

    % 2. Ensure edges are unique and undirected
    % Sort each row so (u,v) is the same as (v,u)
    edgePairs = sort(edgePairs, 2);
    uniqueEdges = unique(edgePairs, 'rows');

    % 3. Calculate Euclidean distances (Weights)
    % Weight = sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)
    v1 = vertices(uniqueEdges(:, 1), :);
    v2 = vertices(uniqueEdges(:, 2), :);
    
    weights = sqrt(sum((v1 - v2).^2, 2));

    % 4. Create the graph object
    G = graph(uniqueEdges(:, 1), uniqueEdges(:, 2), weights);
    
    % Optional: Store vertex coordinates in the graph for later use
    G.Nodes.X = vertices(:, 1);
    G.Nodes.Y = vertices(:, 2);
    G.Nodes.Z = vertices(:, 3);
end