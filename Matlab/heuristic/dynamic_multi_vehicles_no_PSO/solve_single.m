function [served_node,start_time,infeas_nodes,current_time,load,waiting]=solve_single(served_node_all,served_node_v,infeas_nodes_v,start_time_mat_v,current_time_v,load_vec_v,cap,speed,Max_ride,x_coord,y_coord,d_vec,q_vec,e_vec,l_vec)

%%%%%%%%%%%%% Defult values for output%%%%
served_node=[];
start_time=[];
infeas_nodes=[];
load=[];
waiting=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Updating node numbers dynamic part%%%%
n=(size(q_vec,1)-2)/2;
% if dyn_req==1
%     for i_reord=1:size(served_node_v,1)
%         
%         if served_node_v(i_reord,1)>=n
%             served_node_v(i_reord,1)=served_node_v(i_reord,1)+1
%             
%             
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

current_time=current_time_v;
  
  cap_leaving=load_vec_v;
  
   

  
  serving_node_test=served_node_v(end);
 
 
  served_node_coord=[x_coord(serving_node_test+1,1),y_coord(serving_node_test+1,1)]';
  
  v_cap=cap;
 
 infeas_list=infeas_nodes_v;

 waiting_time=[];
 
  served_nodes=served_node_v;
  

 
 
%%%% from here
     
    
    
    
     [travelling_time,nearest_nodes,infeas_nodes,not_good_neghbor]=closest(x_coord,y_coord,e_vec,l_vec,Max_ride,served_node_coord,served_node_all,current_time,speed,d_vec,infeas_list);

  
     i_nearet=size(nearest_nodes,1);
     
    
    
     
     
     
     %%%%%Giving Priority to passangers on board%%
     Pass_onboard=[];
     i_count=1;
     for i_passon=2:size(served_nodes,1)
         if (served_nodes(i_passon,1)<n+1&&ismember(served_nodes(i_passon,1)+n,served_nodes)==0)
             
             Pass_onboard(i_count,1)=served_nodes(i_passon,1);
             Pass_onboard(i_count,2)=i_passon;
             i_count=i_count+1;
             
         end
     end
    
     i_prio=1;
     priority=[];
     count_prio=[];
    
     while (i_prio<=size(Pass_onboard,1))
         
         dist_prio=sqrt((x_coord(Pass_onboard(i_prio,1)+1,1)-served_node_coord(1))^2+(y_coord(Pass_onboard(i_prio,1)+1,1)-served_node_coord(2))^2)/speed;
         dist_destination=sqrt((served_node_coord(1)-x_coord(Pass_onboard(i_prio,1)+n+1,1))^2+(served_node_coord(2)-y_coord(Pass_onboard(i_prio,1)+n+1,1))^2)/speed;
         
         
         time_prio=current_time+dist_prio+d_vec(Pass_onboard(i_prio,1)+1,1);
         wait_prio=max(e_vec(serving_node_test+1,1)-time_prio,0);
        req_time=dist_prio+d_vec(serving_node_test+1,1)+wait_prio;
         Ride_viol=current_time-start_time_mat_v(Pass_onboard(i_prio,2),1)+req_time;
         
         
         time_dest=current_time+dist_destination+req_time;
         wait_dest=max(e_vec(Pass_onboard(i_prio,1)+1+n,1)- time_dest,0);
        Ride_viol_dest=Ride_viol+dist_destination+wait_dest;
        Late_arrival=(current_time+req_time+dist_destination-l_vec(Pass_onboard(i_prio,1)+n+1,1));
       
        
        if serving_node_test<n+1
             dist_dest_dest2=sqrt((x_coord(serving_node_test+n+1)-x_coord(Pass_onboard(i_prio,1)+n+1,1))^2+(y_coord(serving_node_test+n+1)-y_coord(Pass_onboard(i_prio,1)+n+1,1))^2)/speed;
          Late_arrival2=(current_time+req_time+dist_destination+wait_dest+dist_dest_dest2+d_vec(Pass_onboard(i_prio,1)+n+1,1)-l_vec(serving_node_test+n+1,1));
        else
            Late_arrival2=-1;
         end
         
         
         if (Late_arrival>0)||(time_dest>l_vec(Pass_onboard(i_prio,1)+n+1,1))||i_nearet==1||(Late_arrival2>0)%||(Ride_viol> Max_ride)||(Ride_viol_dest>Max_ride)
             priority=[priority;Pass_onboard(i_prio,1)];
             
             count_prio=[count_prio;Pass_onboard(i_prio,2)];
         end
         i_prio=i_prio+1;
     end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if i_nearet>0
         serving_node_test=nearest_nodes(1);
         
         
         i_node=1;
         required_time=d_vec(serving_node_test+1,1)+travelling_time(i_node,1);
         check_prio=0;
    %%%%%%%%Excluding priority%%%
    %priority=[];
    %%%%%%%%%%%%%%%%
         if (size(priority,1)>=1)
             serving_node_test=priority(1,1)+n;
             
             
             travel_time_priority=sqrt((x_coord(served_node_v(end)+1,1)-x_coord(priority(1,1)+n+1,1))^2+(y_coord(served_node_v(end)+1,1)-y_coord(priority(1,1)+n+1,1))^2)/speed;
             check_prio=1;
             
             
             
         else

         
             while ((serving_node_test<n+1&&cap_leaving+q_vec(serving_node_test+1)>v_cap)||((serving_node_test>n)&&ismember(serving_node_test-n,served_nodes)==0)||(serving_node_test<n+1&&ismember(serving_node_test,not_good_neghbor)==1))||((serving_node_test>n)&&ismember(serving_node_test,served_nodes)==1)
                 if i_node<size(nearest_nodes,1)
                 i_node=i_node+1;
                 
                 serving_node_test=nearest_nodes(i_node);
                 
                 
                 
                 required_time=d_vec(serving_node_test+1,1)+travelling_time(i_node,1);
                 else
                    display('No further ')
                    serving_node_test=[];
                    current_time=[];
                     break
                 end
                 
             end
             
%          
         end
         
        
         served_node=serving_node_test;
         
         if check_prio==1&&size(serving_node_test,1)>0
             current_time=current_time+d_vec(served_node_v(end)+1,1)+travel_time_priority;
         elseif size(serving_node_test,1)>0
             current_time=current_time+d_vec(serving_node_test+1,1)+travelling_time(i_node,1);
         end
         
      if size(current_time,1)>0
         if current_time< e_vec(serving_node_test+1,1)
              waiting_time(1,1)=serving_node_test;
               waiting_time(1,2)=e_vec(serving_node_test+1,1)-current_time;
               current_time=current_time+waiting_time(1,2);
             
         end
      else
          waiting_time=[];
      end
         if size(serving_node_test,1)>0
             
             cap_leaving=cap_leaving+q_vec(serving_node_test+1,1);
             
             
             
             
            
             served_node_coord=[x_coord(serving_node_test+1,1),y_coord(serving_node_test+1,1)]';
             start_time=current_time;
             load=cap_leaving;
         end
         
         elseif i_nearet==0
            serving_node_test=[];
            display('No further feasible nearst point')
            i_nearet=0;
            no_neighbor=1;
     
     end
     
    % delay_anim=1;
    %%%if still some passenger onboard%%%
     if i_nearet==0&&size(Pass_onboard,1)>0 
         while size(Pass_onboard,1)>0
       
         served_node=Pass_onboard(1,1)+n;
         
         current_time=current_time+sqrt((served_node_coord(1)-x_coord(Pass_onboard(1,1)+n+1,1))^2+(served_node_coord(2)-y_coord(Pass_onboard(1,1)+n+1,1))^2)/speed+d_vec(served_node_v(end)+1,1);
      
         served_node_coord=[x_coord(Pass_onboard(1,1)+n+1), y_coord(Pass_onboard(1,1)+n+1,1)];
         start_time=current_time;
         load=cap_leaving+q_vec(Pass_onboard(1,1)+n+1,1);
         
         Pass_onboard(1,:)=[];
         end
     end

     waiting=waiting_time;


  
    