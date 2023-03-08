function [remaining,accepted_node,infeasible_node,update_time,leaving_load,waiting_time]=ddarp_one(x_coord,y_coord,e_vec,l_vec,q_vec,v_cap,Max_ride,served_node_coord,served_nodes,served_current_vh,current_time,speed,d_vec,infeas_list,cap_leaving)

           

[travelling_time,nearest_nodes,infeas_nodes,not_good_neghbor]=closest_multi(x_coord,y_coord,e_vec,l_vec,Max_ride,served_node_coord,served_nodes,current_time,speed,d_vec,infeas_list);

n=(size(e_vec,1)-2)/2;
 i_dyn=0;

 infeasible_node=infeas_nodes;
 %%%%Adding infeas node incase it was already visited%%
 [sokhma_node, yeri]=ismember(infeas_nodes,served_current_vh);
 
 if size(sokhma_node,1)>0
    
    for i_s=1:size(sokhma_node,1)
        if yeri(i_s,1)~=0
            node_sokh=served_current_vh(yeri(i_s,1))+n;
           nearest_nodes=[node_sokh;nearest_nodes];
           dist_sok=sqrt((x_coord( node_sokh+1,1)-served_node_coord(1))^2+(y_coord( node_sokh+1,1)-served_node_coord(2))^2)+d_vec(served_current_vh(end)+1,1);
            travelling_time=[dist_sok;travelling_time];
        end
    
    end
     
 end
 i_nearet=size(nearest_nodes,1); 

 %%%%%%%%
 if i_nearet>0
                  serving_node_test=nearest_nodes(1);
                  i_node=1;
                  required_time=d_vec(serving_node_test+1,1)+travelling_time(i_node,1);
                  check_prio=0;
%                   if (size(priority,1)>=1)
%                       serving_node_test=priority(1,1)+n;
%                       
%                       
%                       travel_time_priority=sqrt((x_coord(priority(1,1)+1,1)-x_coord(priority(1,1)+n+1,1))^2+(y_coord(priority(1,1)+1,1)-y_coord(priority(1,1)+n+1,1))^2)/speed;
%                       check_prio=1;
%                       
                      
                      
%                   else
                    
                  i_shombol=0;
                  
                      while ((serving_node_test<n+1&&cap_leaving+q_vec(serving_node_test+1)>v_cap)||((serving_node_test>n)&&ismember(serving_node_test-n,served_current_vh)==0)||(serving_node_test<n+1&&ismember(serving_node_test,not_good_neghbor)==1))
                         
                          if i_node>i_nearet
                             serving_node_test=2*n+1;
                             i_shombol=1;
%                              
                             break
                          
                              
                          else
                          serving_node_test=nearest_nodes(i_node);
                          required_time=d_vec(serving_node_test+1,1)+travelling_time(i_node,1);
                          i_node=i_node+1;
                          end
                          
                      end
                      accepted_node=serving_node_test;
                      if check_prio==1&&size(serving_node_test,1)>0
                          update_time=current_time+d_vec(serving_node_test+1,1)+travel_time_priority;
                      elseif i_shombol==1&&size(serving_node_test,1)>0
                          
                          x_fin=x_coord(served_current_vh(end)+1);
                 y_fin=y_coord(served_current_vh(end)+1);
                 served_node_coord=[x_fin,y_fin];
                 dist_final=sqrt((served_node_coord(1)-x_coord(2*n+2,1))^2+(served_node_coord(2)-y_coord(2*n+2,1))^2)/speed;
                 update_time=current_time+dist_final+d_vec(served_current_vh(end)+1);
                          
                      elseif size(serving_node_test,1)>0
                          update_time=current_time+d_vec(serving_node_test+1,1)+travelling_time(i_node,1);
                      end
                      waiting_time(1,1)=serving_node_test;
                  if  update_time< e_vec(serving_node_test+1,1)
                      
                      waiting_time(1,2)=e_vec(serving_node_test+1,1)-current_time;
                      update_time=update_time+waiting_time(1,2);
                     
                  else
                      waiting_time(1,2)=0;
                  end
                      
                  
                  if size(serving_node_test,1)>0
                      
                      leaving_load=cap_leaving+q_vec(serving_node_test+1,1);
                      
                      
                      
                      
                      
                  end
                  elseif i_nearet==0
                  serving_node_test=[];
                  display('No further feasible nearst point')
                 accepted_node=2*n+1;
                 x_fin=x_coord(served_current_vh(end)+1);
                 y_fin=y_coord(served_current_vh(end)+1);
                 served_node_coord=[x_fin,y_fin];
                 dist_final=sqrt((served_node_coord(1)-x_coord(2*n+2,1))^2+(served_node_coord(2)-y_coord(2*n+2,1))^2)/speed;
                 update_time=current_time+dist_final+d_vec(served_current_vh(end)+1);
                 leaving_load=0;
                 waiting_time=[served_current_vh(end),0];
                 infeasible_node=[]; 
 end
            
 

  remaining=size(nearest_nodes,1);               