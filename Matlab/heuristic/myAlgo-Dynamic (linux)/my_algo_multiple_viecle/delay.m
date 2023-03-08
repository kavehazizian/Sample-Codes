 function F_out=delay(ride_time,waiting_time_mat,Max_ride,l_vec,start_time_mat,served_nodes,n)
 F_vec=inf*ones(size(served_nodes,1),1);
 i_zig=1;

 zig=0;
 for i=1:size(served_nodes,1)
     if size(waiting_time_mat,1)>0
         if served_nodes(i,1)==waiting_time_mat(i_zig,1)&&i_zig<size(waiting_time_mat,1)
             
             zig=zig+waiting_time_mat(i_zig,2);
             i_zig=i_zig+1;
         end
     end
     if served_nodes(i,1)>n&&served_nodes(i,1)~=2*n+1
         i_ride=1;
         while served_nodes(i,1)-n~=ride_time(i_ride,1)
             
             i_ride=i_ride+1;
         end
         P_j=ride_time(i_ride,2);
     else
         P_j=0;
     end
    
     other_term=min(l_vec(served_nodes(i,1)+1,1)-start_time_mat(i),Max_ride-P_j);
     F_vec(i)=zig+max(other_term,0);
 end
 F_out=min(min(F_vec),zig);