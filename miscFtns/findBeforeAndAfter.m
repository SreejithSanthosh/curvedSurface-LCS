function [t1_idX_b,t1_idX_a,wt1,wt2] = findBeforeAndAfter(time,Time_ct)

tol = 10^(-10); %This is so that if the time asked is on the array only have to get one

if (Time_ct-time(end))>tol || -(Time_ct-time(1))>tol
    disp('The time asked for is beyond the range')
end

for i=1:numel(time)-1
    if Time_ct<=time(i+1) && Time_ct>=time(i)
        t1_idX_b = i;
        t1_idX_a = i+1;
        
        dist2_b = Time_ct-time(i);
        dist2_a = time(i+1)-Time_ct;
        
        if dist2_b<tol
            wt2 = 0;
            wt1 = 1;
        elseif dist2_a<tol
            wt1 = 0;
            wt2 = 1;
        else
            wt1 = (time(i+1)-Time_ct)/(time(i+1)-time(i));
            wt2 = (Time_ct-time(i))/(time(i+1)-time(i));
        end 
    end 
end 

% This is for the cases where time slipped slightly off min and max 
if abs(Time_ct-time(1))<tol && time(1)-Time_ct>0
    wt1 = 1;
    wt2 = 0;
    t1_idX_b = 1;
    t1_idX_a = 1+1;
end 

if abs(Time_ct-time(end))<tol && time(end)-Time_ct<0
    wt1 = 0;
    wt2 = 1;
    t1_idX_b = numel(time)-1;
    t1_idX_a = numel(time);
end 


if wt1<0 ||wt2<0
    disp('The weights went negative')
end
    
    
end 