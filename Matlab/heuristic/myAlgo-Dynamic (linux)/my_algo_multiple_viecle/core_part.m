
function [starting_time,Total_dist_modified]=core_part(alpha,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec,n_s,m,max_route, Max_ride,v_cap,v_routing_time,v_speed,initial_time)
global_time=0;
 pool_size=2*n_s+2;
 time_res=5;
    for i=1:m
    served_nodes{i}=[0];
    start_time_mat{i}=initial_time(i);
    infeas_nodes_mat{i}=[];
     load_vec{i}=[0];
     current_time_mat{i}=initial_time(i);
     waiting_time_mat{i}=zeros(1,2);
    end
    
 served_node_all=[0];
 
 

 
 %%%%%%%%Initiating some variables%%%%

 
  n_d=0;
 i_dyn=1;
                                                                                                                                                                    
 n=n_s+n_d;
 n_add=n_s;
time_windows_viol=0;
total_ride_viol=0;
route_duration=0;
%total_waiting=0;
%route_viol=0;
 gamma=4;
 %%%%%%%%%%%%%%%%%%%%%%%%


while (global_time<max_route&&pool_size>0)
     
     for i_v=1:m
                           %%%% current_time_mat{i_v}(1) allows to start serving even after time zero%%%%   
         if current_time_mat{i_v}(end)<=global_time+current_time_mat{i_v}(1)
             cap=v_cap(i_v);
             rout_time=v_routing_time(i_v); 
             speed=v_speed(i_v);
          
             [served_node,start_time,infeas_nodes,current_time,load,waiting]=solve_single(served_node_all,served_nodes{i_v},infeas_nodes_mat{i_v},start_time_mat{i_v},current_time_mat{i_v}(end),load_vec{i_v}(end),cap,speed,Max_ride,x_coord,y_coord,d_vec,q_vec,e_vec,l_vec,alpha);
            
             
             
             served_nodes{i_v}=[served_nodes{i_v};served_node];
             start_time_mat{i_v}=[start_time_mat{i_v};start_time];
             infeas_nodes_mat{i_v}=[infeas_nodes_mat{i_v};infeas_nodes];
             current_time_mat{i_v}=[current_time_mat{i_v};current_time];
             load_vec{i_v}=[load_vec{i_v};load];
             
             waiting_time_mat{i_v}=[waiting_time_mat{i_v};waiting];
             if served_node>0
                 served_node_all=[served_node_all;served_node];
             end
             
          
             
          
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             
             
             
             
             
             
             
             
             
             
             
         end
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         
         
         

         
     end

     
 
        
          
          
          
          
          %%%%%%%%%%%%%Updating the number of served nodes and infeasible nodes%%%
          
          
%           for i_reordppp=2:size(served_node_all,1)
%               if served_node_all(i_reordppp,1)>=n_add
%                   
%                   served_node_all(i_reordppp,1)=served_node_all(i_reordppp,1)+1;
%                   
%               end
%           end
%           
%           i_dyn=i_dyn+1;
%           
%           
%       end
      
      
      
      
      
      
      global_time=global_time+time_res;
      pool_size=2*n_s+1+2*n_d-size(served_node_all,1);
      
 end
 
 
 %%%%%%%%%%%%%%%%%%%last line to deopt%%%%%%%%%%%%%%%%%%%

 for i_l=1:m
     
     %%%%%%%%%%%%%%%%%%%%Allowing the last picked-up customer to reach iits
     %%%%%%%%%%%%%%%%%%%%destination%%%%
     i_rem_co=size(served_nodes{i_l},1);
     still_onboard=[];
     for i_pass=i_rem_co:-1:2
         if ismember(served_nodes{i_l}(i_pass)+n,served_nodes{i_l})==0&&served_nodes{i_l}(i_pass)<n+1
             still_onboard=[still_onboard;served_nodes{i_l}(i_pass)+n];
             
             
         end
     end
     
     
     if still_onboard~=0
         for i_si=size(still_onboard,1):-1:1
             
             served_nodes{i_l}=[served_nodes{i_l};still_onboard(i_si)];
             
             x_fin=x_coord(served_nodes{i_l}(end-1)+1);
             y_fin=y_coord(served_nodes{i_l}(end-1)+1);
             served_node_coord=[x_fin,y_fin];
             dist_final=sqrt((served_node_coord(1)-x_coord(served_nodes{i_l}(end-2),1))^2+(served_node_coord(2)-y_coord(served_nodes{i_l}(end-2),1))^2)/v_speed(i_l,1);
             load_vec{i_l}=[load_vec{i_l};load_vec{i_l}(end)+q_vec(served_nodes{i_l}(end-1)+1,1)];
             start_time_mat{i_l}(end+1)=start_time_mat{i_l}(end)+dist_final+d_vec(served_nodes{i_l}(end-1)+1);
             
             current_time=start_time_mat{i_l}(end);

         end
         
     end
     
     %%%%%%%%%%%%%%%%%%%%%
     
     if served_nodes{i_l}(end)~=2*n+1;
         served_nodes{i_l}=[served_nodes{i_l};2*n+1];
         x_fin=x_coord(served_nodes{i_l}(end-1)+1);
         y_fin=y_coord(served_nodes{i_l}(end-1)+1);
         served_node_coord=[x_fin,y_fin];
         dist_final=sqrt((served_node_coord(1)-x_coord(2*n+2,1))^2+(served_node_coord(2)-y_coord(2*n+2,1))^2)/v_speed(i_l,1);
         load_vec{i_l}=[load_vec{i_l};0];
         start_time_mat{i_l}(end+1)=start_time_mat{i_l}(end)+dist_final+d_vec(served_nodes{i_l}(end-1)+1);
         
         current_time=start_time_mat{i_l}(end);
         

     else
         while served_nodes{i_l}(end)==served_nodes{i_l}(end-1)
             served_nodes{i_l}(end)=[];
             start_time_mat{i_l}(end)=[];
             load_vec{i_l}(end)=[];
         end
     end
     %%%%Calculating the ride  time and excess ride time%%%%
     [ride_time{i_l},ride_viol(i_l)]=ride_time_cal(served_nodes{i_l},start_time_mat{i_l},d_vec,Max_ride);
     total_ride_viol=total_ride_viol+ride_viol(i_l);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%Calculating the violence of constraints%%%%
 time_windows_viol=time_windows_viol+time_viol_cal(served_nodes{i_l},start_time_mat{i_l},l_vec);
 %%%%%max_route_violence%%%
% route_viol=route_viol+max(start_time_mat{i_l}(end)-v_routing_time(i_l),0);
 route_duration=route_duration+(start_time_mat{i_l}(end)-start_time_mat{i_l}(1));
 %%%%%%%%%Calculating waiting time sum%%
 %total_waiting=total_waiting+sum(waiting_time_mat{i_l}(:,2));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     
     
 end
 
 
 tot_dist=zeros(m,1);
 starting_time=zeros(m,1);
 for i=1:m
     tot_dist(i)=distance_travelled(x_coord,y_coord,served_nodes{i});
   F_0=delay(ride_time{i},waiting_time_mat{i},Max_ride,l_vec,start_time_mat{i},served_nodes{i},n);
     starting_time(i)=e_vec(1,1)+F_0;
 end
 Total_distance_traveled=sum(tot_dist);
 beta=2*n/(size(served_node_all,1)-1);
 Total_dist_modified=-(1+beta^gamma*(beta-1))* Total_distance_traveled-gamma*(time_windows_viol+total_ride_viol);%-route_duration;%)+route_viol);%-total_waiting;

 %%%%%%% Loop for solving darp ends here%%%%%%%%%%%