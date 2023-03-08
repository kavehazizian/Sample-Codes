
function [starting_time,Total_dist_modified]=core_part2(alpha,x_coord,y_coord,d_vec, q_vec,e_vec,...
    l_vec,m,max_route, Max_ride,v_cap,v_routing_time,v_speed,req_time,...
    initial_depots_spec,served_nodesd,start_time_matd,load_vecd,waiting_time_matd)

n_s=(size(e_vec,1)-1)/2;

global_time=req_time;
 pool_size=2*n_s+2;
 time_res=5;
 served_node_all=[0];
    for i=1:m
    served_nodes{i}=served_nodesd{i};
    start_time_mat{i}=start_time_matd{i};
    infeas_nodes_mat{i}=[];
     load_vec{i}=load_vecd{i};
     current_time_mat{i}=sort(start_time_matd{i});
     waiting_time_mat{i}=waiting_time_matd{i};
     served_node_all=[served_node_all;served_nodesd{i}(2:end)];
    end
    

 
 %%%%%%%%%%%Maintatinig initia inputs%%%%
 x_coordp=x_coord;
 y_coordp=y_coord;
 d_vecp=d_vec;
 q_vecp=q_vec;
 e_vecp=e_vec;
 l_vecp=l_vec;
 
 %%%%%%%%%

 
 %%%%%%%%Initiating some variables%%%%

                                                                                                                                                                     
 n=n_s;
 
time_windows_viol=0;
total_ride_viol=0;
total_waiting=0;
 gamma=4;
 %%%%%%%%%%%%%%%%%%%%%%%%


while (global_time<max_route&&pool_size>0)
     
     for i_v=1:m
        %%%%%%%Adding inital depots%%%%%
         x_coord=[initial_depots_spec(i_v,1);x_coordp];
         y_coord=[initial_depots_spec(i_v,2);y_coordp];
         d_vec=[initial_depots_spec(i_v,3);d_vecp];
         q_vec=[initial_depots_spec(i_v,4);q_vecp];
         e_vec=[initial_depots_spec(i_v,5);e_vecp];
         l_vec=[initial_depots_spec(i_v,6);l_vecp];
         
         
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                           %%%% current_time_mat{i_v}(1) allows to start serving even after time zero%%%%   
         if current_time_mat{i_v}(end)<=global_time%+current_time_mat{i_v}(1)
             cap=v_cap(i_v);
             rout_time=v_routing_time(i_v); 
             speed=v_speed(i_v);
          
             [served_node,start_time,infeas_nodes,current_time,load,waiting]=solve_single(served_node_all,served_nodes{i_v},infeas_nodes_mat{i_v},start_time_mat{i_v},...
                 current_time_mat{i_v}(end),load_vec{i_v}(end),cap,speed,Max_ride,x_coord,y_coord,d_vec,q_vec,e_vec,l_vec,alpha);
            
             
             
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

     
 
      
      
      
      
      
      
      global_time=global_time+time_res;
      pool_size=2*n_s+1-size(served_node_all,1);
      
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

             x_fin=x_coord(served_nodes{i_l}(end)+1);
             y_fin=y_coord(served_nodes{i_l}(end)+1);
             served_node_coord=[x_fin,y_fin];
             
             dist_final=sqrt((served_node_coord(1)-x_coord(served_nodes{i_l}(end-1)+1,1))^2+(served_node_coord(2)-y_coord(served_nodes{i_l}(end-1)+1,1))^2)/v_speed(i_l,1);
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
 
 %%%%%%%%%Calculating waiting time sum%%
 if size (waiting_time_mat{i_l},1)>0
 total_waiting=total_waiting+sum(waiting_time_mat{i_l}(:,2));
 end
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
 Total_dist_modified=-(1+beta^gamma*(beta-1))* Total_distance_traveled-gamma*(time_windows_viol+total_ride_viol);%-total_waiting;

 %%%%%%% Loop for solving darp ends here%%%%%%%%%%%