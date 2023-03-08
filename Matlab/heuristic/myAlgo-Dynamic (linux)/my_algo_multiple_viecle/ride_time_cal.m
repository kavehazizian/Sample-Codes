function [ride_time,ride_viol]=ride_time_cal(served_nodes,start_time_mat,d_vec,Max_ride)
i_r=0;
n=(size(d_vec,1)-2)/2;
ride_time=[];

for i=2:size(served_nodes,1)-1
    if served_nodes(i,1)<n+1
        ip=i+1;
       
        while ((served_nodes(i,1)+n)~=served_nodes(ip,1))
            
            ip=ip+1;
        end
        i_r=i_r+1;
      ride_time(i_r,1)=served_nodes(i,1); 
      ride_time(i_r,2)=start_time_mat(ip,1)-start_time_mat(i,1)-d_vec(served_nodes(i,1)+1,1);
        
    end
end
no_served=size(ride_time,1);
viol_vec=zeros(no_served,1);
if no_served>0
    for i_v=1:no_served
        diff=ride_time(i_v,2)-Max_ride;
        if diff>0
            viol_vec(i_v,1)=diff;
        end
    end
    ride_viol=sum(viol_vec(:,1));
else
    
    ride_viol=0;
end