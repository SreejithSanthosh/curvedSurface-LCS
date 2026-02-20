function [sv_list,vec0_list,vecf_list] = compute_deform_euclidean(r0,rf,delta)
% r0: (Nq x d) double array
% rf: (Nq x d) double array
% delta : scalar double : spatial smoothing
% d is the dimension

% Output 
% sv_list : Nq x d: sv_list(:,1)< sv_list(:,2)<..
% vec0_list : Nq x d x d: vec0_list(i,:,1) - vec0_list(i,:,2) -..
% vecf_list : Similar to vec0_list
 
dim = size(r0,2); % Is it 2D or 3D (d = 2/3) 
Nq = size(r0,1); % Number of points being advected

% Intialize the variables
sv_list = nan(Nq, dim); % Singular value of DF
vec0_list = nan(Nq, dim,dim); % Eig. vector of DF at x0
vecf_list = nan(Nq, dim,dim); % Eig. vector of DF at xf

h = waitbar(0);
for i = 1:Nq
    % Computes the nearest nodes within this distance
    neighb = vecnorm(r0-r0(i,:),2,2)<delta; neighb(i) =0;
    
    % Consitruct the matrix 
    X = r0(neighb,:)-r0(i,:); Y = rf(neighb,:)-rf(i,:);

    warning('off','MATLAB:rankDeficientMatrix')
    DF = Y'/X'; % (3 x 3) double  
    warning('on','MATLAB:rankDeficientMatrix')
        

    [V,D] = eig(DF'*DF);
    [D,I] = sort(diag(D)); % Sort in ascending order  
    
    V0 = V(:,I); % Makes the matrix is V0: (d x d) :eigVec at x0
    Vf = DF*V0; %  Makes the matrix is Vf: (d x d) :eigVec at xf

    sv_list(i,:) = sqrt(D);    
    vec0_list(i,:,:) = V0;
    vecf_list(i,:,:) = Vf;

    % Update the waitbar
    waitbar(i/Nq, h);
end 

delete(h); % Close the bar when finished

end 
