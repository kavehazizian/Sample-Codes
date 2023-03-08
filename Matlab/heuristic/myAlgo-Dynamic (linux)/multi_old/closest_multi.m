function [travelling_time,nearest_nodes,infeas_nodes,not_good_neghbor]=closest_multi(x_coord,y_coord,e_vec,l_vec,Max_ride,served_node_coord,served_nodes,current_time,speed,d_vec,infeas_list)
n=(size(e_vec,1)-2)/2;

node_no=served_nodes(end);


x_current=served_node_coord(1);
y_current=served_node_coord(2);
if node_no<n+1
   x_dest=x_coord(node_no+n+1,1); 
   y_dest=y_coord(node_no+n+1,1);


end



distance=[];
travelling_time_or=zeros(2*n+1-size(served_nodes,1),1);
travelling_time=zeros(2*n+1-size(served_nodes,1),1);
infeas_nodes=[];
not_good_neghbor=[];

alpha=1;
beta=1;
gamma=1;
i_node=1;
i_inf=1;
for i=1:2*(n)
    if (ismember(i,served_nodes)==0)&&(ismember(i,infeas_list)==0)
        dist=sqrt((x_coord(i+1,1)-x_current)^2+(y_coord(i+1,1)-y_current)^2);
        
       travelling_time_or(i_node,1)=dist/speed;
       arrival_time=current_time+travelling_time_or(i_node,1)+d_vec(node_no+1,1);
       waiting_time=e_vec(i+1,1)-arrival_time;
       
       
       
       
       if i<n+1
           %%%From current node to candidate node%%%
          dist_dest=sqrt((x_coord(i+1,1)-x_coord(i+n+1,1))^2+(y_coord(i+1,1)-y_coord(i+n+1,1))^2); 
          travel_time_dest=dist_dest/speed;
          
          arrival_dest=arrival_time+max(waiting_time,0)+d_vec(i+1,1)+travel_time_dest;
          
          feas_check=l_vec(i+n+1,1)-arrival_dest;
          
          
          if node_no<n+1
       %%%%%from candidate node to destination of the current node if it is
       %%%%%a pickup node
              dest_dist_c=sqrt((x_coord(i+1,1)-x_dest)^2+(y_coord(i+1,1)-y_dest)^2);
       travel_tim_dest=dest_dist_c/speed;
       arrival_dest_current=arrival_time+travel_tim_dest+d_vec(i+1,1)+max(waiting_time,0);
       waiting_in_des=max(e_vec(node_no+1)-arrival_dest_current,0);
      
       %%%%From destination of current node to destination of candidate
       %%%%node if both r pick up nodes
       dest_to_des=sqrt((x_coord(i+n+1,1)-x_dest)^2+(y_coord(i+n+1,1)-y_dest)^2)/speed;
        arrival_f_des_t_des=arrival_dest_current+ d_vec(node_no+n+1,1)+waiting_in_des+dest_to_des;
     feas_check_des=l_vec(node_no+n+1,1)-arrival_dest_current;
      feas_check_des_des=l_vec(i+n+1,1)-arrival_f_des_t_des;
     if feas_check_des<0||feas_check_des_des<0
        not_good_neghbor=[not_good_neghbor;i];
         
     end
     
          end
      
       end
       
       if l_vec(i+1,1)>=arrival_time&&i>n
           distance(i_node,1)=i;
        distance(i_node,2)=alpha*dist+beta*max(waiting_time,0);
        i_node=i_node+1;
       elseif l_vec(i+1,1)>=arrival_time&&i<n+1&&feas_check>=0
          distance(i_node,1)=i;
        distance(i_node,2)=alpha*dist+beta*max(waiting_time,0)+gamma*max(e_vec(i+n+1)-arrival_dest,0);
        i_node=i_node+1; 
       
           
       else
        
        infeas_nodes(i_inf,1)=i;
        
        i_inf=i_inf+1;
       end
        
    end
end
if size(infeas_nodes,1)>0
for i_ch=1:size(infeas_nodes,1)
    if infeas_nodes(i_ch,1)<n+1&&(ismember(infeas_nodes(i_ch,1)+n,infeas_nodes)==0)
        infeas_nodes=[infeas_nodes;infeas_nodes(i_ch,1)+n];
    elseif infeas_nodes(i_ch,1)>=n+1&&(ismember(infeas_nodes(i_ch,1)-n,infeas_nodes)==0)
        infeas_nodes=[infeas_nodes;infeas_nodes(i_ch,1)-n];
    end
    
end
end






if size(distance,1)>0
[int_vec,~,ind_dist]=intersect(infeas_nodes,distance(:,1));
if size(int_vec,2)>0
for i_del=1:size(ind_dist,1)
distance(ind_dist(i_del),2)=-1;
end
end

[dist_vec,nearest_nodes_or]=sort(distance(:,2));
 
i_d=1;
while (size(nearest_nodes_or,1)>1&&dist_vec(i_d)==-1&&i_d<=size(dist_vec,1))
 
    nearest_nodes_or(1)=[];
 
    
        i_d=i_d+1;
    
end

if size(nearest_nodes_or,1)==1&&dist_vec(1)==-1
     nearest_nodes_or(1)=[];
end


if size(nearest_nodes_or,2)>0
for j=1:size(nearest_nodes_or,1)
   
    nearest_nodes(j,1)=distance(nearest_nodes_or(j,1),1);
    
    travelling_time(j,1)=travelling_time_or(nearest_nodes_or(j,1),1);
end
else
 nearest_nodes=[];
     travelling_time=[];   
end
else
    nearest_nodes=[];
     travelling_time=[];


end

 
