% This is constructed on 1/3/26
function rt_arr = advect_for(mesh_F,mesh_time,mesh_v,mesh_r,r0)

% Fill this out to check 
arguments
   mesh_F
   mesh_time
   mesh_v
   mesh_r
   r0 = mesh_r(1,:) % (Cell 1x3 )This sets the intial condition to nodes at t0;if noinput
end 
% Assuming uniform time and time-stepping through all
Nt = numel(mesh_time); dt = mesh_time(2)-mesh_time(1); 

% Set initial condition 
rt_arr = cell(Nt,3);
rt_arr{1,1} = r0{:,1};
rt_arr{1,2} = r0{:,2};
rt_arr{1,3} = r0{:,3};


for ct = 1:Nt-1
    rp = cell2mat(rt_arr(ct,:));
    
    mesh_F_ct = mesh_F{ct};
    mesh_r_ct = cell2mat(mesh_r(ct,:)); % Since it is multi-dim
    mesh_v_ct = cell2mat(mesh_v(ct,:)); % Since it is multi-dim

    mesh_F_nt = mesh_F{ct+1};
    mesh_r_nt = cell2mat(mesh_r(ct+1,:));
    mesh_v_nt = cell2mat(mesh_v(ct+1,:));

    % The first step in the Heuns method 
    v_ct = interp_on_mesh(mesh_F_ct,mesh_r_ct,mesh_v_ct,rp);
    rpp = rp+dt*v_ct; 

    % The second step in the Heuns method 
    % Dont need to project inbetween as interp does it 
    v_nt = interp_on_mesh(mesh_F_nt,mesh_r_nt,mesh_v_nt,rpp);
    r_nt = rp+(dt/2)*(v_ct+v_nt);   

    % Store the data 
    rt_arr{ct+1,1} = r_nt(:,1);
    rt_arr{ct+1,2} = r_nt(:,2);
    rt_arr{ct+1,3} = r_nt(:,3);
end 

end 