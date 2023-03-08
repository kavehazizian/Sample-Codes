function ride_time=ride_time_cal(served_nodes,start_time_mat,n,d_vec)
i_r=0;
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