function [xt_Advect,yt_Advect,zt_Advect] = AdvctBck(cpu_num,tf,tSave_Advect,tArr_Advect,xq,yq,zq,time,v,x,y,z,TrianT)

dt = tArr_Advect(2) - tArr_Advect(1); Nt = numel(tArr_Advect);

%Project the query points to the mesh 
Time_ct = tf-tArr_Advect(1);
[t1_idX_b,t1_idX_a,wt1,wt2] = findBeforeAndAfter(time,Time_ct);

if wt1>wt2
    t_idX_closest = t1_idX_b;
else
    t_idX_closest = t1_idX_a;
end

xMesh = squeeze(x{t_idX_closest}); yMesh = squeeze(y{t_idX_closest}); zMesh = squeeze(z{t_idX_closest});

inputs = struct(); inputs.nodes = [squeeze(xMesh), squeeze(yMesh),squeeze(zMesh)]; %Faces
inputs.faces = TrianT{t_idX_closest}; points = [xq,yq,zq];

surface_points = fastPoint2TriMeshModif(inputs,points);
xq = squeeze(surface_points(:,1)); yq = squeeze(surface_points(:,2)); zq = squeeze(surface_points(:,3));
    
%Intialize the variables
xt_Advect = zeros(numel(tSave_Advect),numel(xq)); yt_Advect = xt_Advect; zt_Advect= xt_Advect;

FV=struct('faces',{},'vertices',{});

%% Split the initial condition 
N_split = cpu_num;
% d = numel(xq)/N_split;
xSplit = splitArray(xq,N_split); ySplit = splitArray(yq,N_split); zSplit = splitArray(zq,N_split);

%Store the intermediate result in a cell array
xTrajec_cell = cell(N_split,1); yTrajec_cell = cell(N_split,1); zTrajec_cell = cell(N_split,1);

% addAttachedFiles(gcp,'findBeforeAndAfter.m');
parfor k = 1:N_split
        %Build storage for each k
        numPts_k = size(xSplit{k},1);
        xTrajec_k = zeros([numel(tSave_Advect),numPts_k]); yTrajec_k = xTrajec_k;
        zTrajec_k = xTrajec_k;
        
        % Initialize the first time step
        cp = 1;
        % xp = xq(1+d*(k-1):d*(k)); yp = yq(1+d*(k-1):d*(k)); zp = zq(1+d*(k-1):d*(k)); 
        xp = xSplit{k}; yp = ySplit{k}; zp = zSplit{k};
        xTrajec_k(cp,:) = xp; yTrajec_k(cp,:) = yp; zTrajec_k(cp,:) = zp; 

        % progressbar
        for ct = 1:Nt-1
            Time_ct = tf-tArr_Advect(ct);

            %Need to find the velocity for t(ct),t(ct)+0.5*dt and t(ct)+dt
            % In the function findBeforeAndAfter if the time lies exactly at a time
            % step we send back 0 for the other Output
            [t1_idX_b,t1_idX_a,wt1_k1,wt2_k1] = findBeforeAndAfter(time,Time_ct);
            [t2_idX_b,t2_idX_a,wt1_k2,wt2_k2] = findBeforeAndAfter(time,Time_ct-dt/2);
            t3_idX_b = t2_idX_b; t3_idX_a = t2_idX_a; wt1_k3 = wt1_k2; wt2_k3= wt2_k2;
            [t4_idX_b,t4_idX_a,wt1_k4,wt2_k4] = findBeforeAndAfter(time,Time_ct-dt);
            
            %% First step RK4 
            xp_k1 = xp; yp_k1 = yp; zp_k1 = zp;
            v1_b = 0*xp; v2_b = 0*xp; v3_b = 0*xp; 

            if wt1_k1>0
                % xData = squeeze(x(:,t1_idX_b)); yData = squeeze(y(:,t1_idX_b)); zData = squeeze(z(:,t1_idX_b));
                % v1Data = squeeze(v(:,1,t1_idX_b)); v2Data = squeeze(v(:,2,t1_idX_b)); v3Data = squeeze(v(:,3,t1_idX_b));
                xData = squeeze(x{t1_idX_b}); yData = squeeze(y{t1_idX_b}); zData = squeeze(z{t1_idX_b});
                v1Data = squeeze(v{1,t1_idX_b}); v2Data = squeeze(v{2,t1_idX_b}); v3Data = squeeze(v{3,t1_idX_b});
                
                [v1_b,v2_b,v3_b] = interpOnMesh(TrianT{t1_idX_b},xData,yData,zData,v1Data,v2Data,v3Data,xp_k1,yp_k1,zp_k1);
            end 

            v1_a = 0*xp; v2_a = 0*xp; v3_a = 0*xp; 

            if wt2_k1>0

                xData = squeeze(x{t1_idX_a}); yData = squeeze(y{t1_idX_a}); zData = squeeze(z{t1_idX_a});
                v1Data = squeeze(v{1,t1_idX_a}); v2Data = squeeze(v{2,t1_idX_a}); v3Data = squeeze(v{3,t1_idX_a});

                [v1_a,v2_a,v3_a] = interpOnMesh(TrianT{t1_idX_a},xData,yData,zData,v1Data,v2Data,v3Data,xp_k1,yp_k1,zp_k1);
            end 

            v1_k1 = wt1_k1*v1_b+wt2_k1*v1_a; v2_k1 = wt1_k1*v2_b+wt2_k1*v2_a; v3_k1 = wt1_k1*v3_b+wt2_k1*v3_a; % The first velocity component

            %% Second step RK4 
            xp_k2 = xp+(dt/2)*v1_k1; yp_k2 = yp+(dt/2)*v2_k1; zp_k2 = zp+(dt/2)*v3_k1;
            v1_b = 0*xp; v2_b = 0*xp; v3_b = 0*xp; 

            if wt1_k2>0
                xData = squeeze(x{t2_idX_b}); yData = squeeze(y{t2_idX_b}); zData = squeeze(z{t2_idX_b});
                v1Data = squeeze(v{1,t2_idX_b}); v2Data = squeeze(v{2,t2_idX_b}); v3Data = squeeze(v{3,t2_idX_b});

                [v1_b,v2_b,v3_b] = interpOnMesh(TrianT{t2_idX_b},xData,yData,zData,v1Data,v2Data,v3Data,xp_k2,yp_k2,zp_k2);
            end 

            v1_a = 0*xp; v2_a = 0*xp; v3_a = 0*xp; 

            if wt2_k2>0

                xData = squeeze(x{t2_idX_a}); yData = squeeze(y{t2_idX_a}); zData = squeeze(z{t2_idX_a});
                v1Data = squeeze(v{1,t2_idX_a}); v2Data = squeeze(v{2,t2_idX_a}); v3Data = squeeze(v{3,t2_idX_a});

                [v1_a,v2_a,v3_a] = interpOnMesh(TrianT{t2_idX_a},xData,yData,zData,v1Data,v2Data,v3Data,xp_k2,yp_k2,zp_k2);
            end 

            v1_k2 = wt1_k2*v1_b+wt2_k2*v1_a; v2_k2 = wt1_k2*v2_b+wt2_k2*v2_a; v3_k2 = wt1_k2*v3_b+wt2_k2*v3_a; % The second velocity component

            %% Third step RK4 
            xp_k3 = xp+(dt/2)*v1_k2; yp_k3 = yp+(dt/2)*v2_k2; zp_k3 = zp+(dt/2)*v3_k2;
            v1_b = 0*xp; v2_b = 0*xp; v3_b = 0*xp; 

            if wt1_k3>0

                xData = squeeze(x{t3_idX_b}); yData = squeeze(y{t3_idX_b}); zData = squeeze(z{t3_idX_b});
                v1Data = squeeze(v{1,t3_idX_b}); v2Data = squeeze(v{2,t3_idX_b}); v3Data = squeeze(v{3,t3_idX_b});

                [v1_b,v2_b,v3_b] = interpOnMesh(TrianT{t3_idX_b},xData,yData,zData,v1Data,v2Data,v3Data,xp_k3,yp_k3,zp_k3);
            end 

            v1_a = 0*xp; v2_a = 0*xp; v3_a = 0*xp; 

            if wt2_k3>0
                xData = squeeze(x{t3_idX_a}); yData = squeeze(y{t3_idX_a}); zData = squeeze(z{t3_idX_a});
                v1Data = squeeze(v{1,t3_idX_a}); v2Data = squeeze(v{2,t3_idX_a}); v3Data = squeeze(v{3,t3_idX_a});

                [v1_a,v2_a,v3_a] = interpOnMesh(TrianT{t3_idX_a},xData,yData,zData,v1Data,v2Data,v3Data,xp_k3,yp_k3,zp_k3);
            end 

            v1_k3 = wt1_k2*v1_b+wt2_k2*v1_a; v2_k3 = wt1_k2*v2_b+wt2_k2*v2_a; v3_k3 = wt1_k2*v3_b+wt2_k2*v3_a; % The third velocity component
            
            %% Fourth step RK4 
            xp_k4 = xp+(dt)*v1_k3; yp_k4 = yp+(dt)*v2_k3; zp_k4 = zp+(dt)*v3_k3;
            v1_b = 0*xp; v2_b = 0*xp; v3_b = 0*xp; 

            if wt1_k4>0
                xData = squeeze(x{t4_idX_b}); yData = squeeze(y{t4_idX_b}); zData = squeeze(z{t4_idX_b});
                v1Data = squeeze(v{1,t4_idX_b}); v2Data = squeeze(v{2,t4_idX_b}); v3Data = squeeze(v{3,t4_idX_b});

                [v1_b,v2_b,v3_b] = interpOnMesh(TrianT{t4_idX_b},xData,yData,zData,v1Data,v2Data,v3Data,xp_k4,yp_k4,zp_k4);
            end 

            v1_a = 0*xp; v2_a = 0*xp; v3_a = 0*xp; 

            if wt2_k4>0
                xData = squeeze(x{t4_idX_a}); yData = squeeze(y{t4_idX_a}); zData = squeeze(z{t4_idX_a});
                v1Data = squeeze(v{1,t4_idX_a}); v2Data = squeeze(v{2,t4_idX_a}); v3Data = squeeze(v{3,t4_idX_a});

                [v1_a,v2_a,v3_a] = interpOnMesh(TrianT{t4_idX_a},xData,yData,zData,v1Data,v2Data,v3Data,xp_k4,yp_k4,zp_k4);
            end 

            v1_k4 = wt1_k2*v1_b+wt2_k2*v1_a; v2_k4 = wt1_k2*v2_b+wt2_k2*v2_a; v3_k4 = wt1_k2*v3_b+wt2_k2*v3_a; % The fourth velocity component
            %% Put all the components together and update the positions 

            xct = xp - (dt/6)*(v1_k1+2*v1_k2+2*v1_k3+v1_k4);
            yct = yp - (dt/6)*(v2_k1+2*v2_k2+2*v2_k3+v2_k4);
            zct = zp - (dt/6)*(v3_k1+2*v3_k2+2*v3_k3+v3_k4);

            % Need to project it down to the mesh closest to time_ct+dt
            if wt1_k4 > wt2_k4
                t_idX_closest = t4_idX_b;
            else 
                t_idX_closest = t4_idX_a;
            end 
            
            % xMesh = squeeze(x(:,t_idX_closest)); yMesh = squeeze(y(:,t_idX_closest)); zMesh = squeeze(z(:,t_idX_closest));
            xMesh = squeeze(x{t_idX_closest}); yMesh = squeeze(y{t_idX_closest}); zMesh = squeeze(z{t_idX_closest});
            vertices = [squeeze(xMesh), squeeze(yMesh),squeeze(zMesh)]; %Faces

            % faces = convhulln([squeeze(xMesh),squeeze(yMesh),squeeze(zMesh)]);
            faces = TrianT{t_idX_closest};

            points = [xct,yct,zct];

            surface_points = point2trimesh_mycompute(vertices,faces,points);
            xct = squeeze(surface_points(:,1)); yct = squeeze(surface_points(:,2)); zct = squeeze(surface_points(:,3));

            %% Decide whether to save the data or not 
            if ismember(tArr_Advect(ct+1),tSave_Advect)
                idX_Save = find(tSave_Advect==tArr_Advect(ct+1));
                xTrajec_k(idX_Save,:) = xct; yTrajec_k(idX_Save,:) = yct; zTrajec_k(idX_Save,:) = zct; 
            end

            %% Set the x_ct variables as the xp variables 
            xp = xct; yp = yct; zp = zct;

            % progressbar(ct/Nt)
        end
        xTrajec_cell{k} = xTrajec_k;  yTrajec_cell{k} = yTrajec_k;  zTrajec_cell{k} = zTrajec_k;
end 

%% Take everything from the cell and put it into the array 

% for k= 1:N_split
%     xt_Advect(:,1+d*(k-1):d*(k)) = xTrajec_cell{k}; yt_Advect(:,1+d*(k-1):d*(k)) = yTrajec_cell{k};
%     zt_Advect(:,1+d*(k-1):d*(k)) = zTrajec_cell{k}; 
% end 

count = 0;
for k = 1:N_split
    numPts_k = size(xSplit{k},1);
    xt_Advect(:,count+1:count+numPts_k) = xTrajec_cell{k}; yt_Advect(:,count+1:count+numPts_k) = yTrajec_cell{k};
    zt_Advect(:,count+1:count+numPts_k) = zTrajec_cell{k}; 
    count = count+numPts_k;
end 
end 

function surface_points = point2trimesh_mycompute(vertices,faces,points)

FV.faces = faces; FV.nodes = vertices;
[surface_points,~,~] = fastPoint2TriMeshModif(FV,points,0);
end