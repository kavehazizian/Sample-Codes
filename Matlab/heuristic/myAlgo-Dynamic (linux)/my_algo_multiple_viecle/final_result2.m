
function [ Total_distance_traveled,Total_dist_modified,served_nodes,start_time_mat,infeas_nodes_mat,load_vec,waiting_time_mat, ride_time]=final_result2(alpha,...
    x_coord,y_coord,d_vec, q_vec,e_vec,l_vec,m,max_route, Max_ride,v_cap,v_routing_time,v_speed,req_time,...
    initial_depots_spec,served_nodesd,start_time_matd,load_vecd,waiting_time_matd)
n_s=(size(e_vec,1)-1)/2;
n=n_s;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  k_mov=ones(m,1);
 
time_windows_viol=0;
total_ride_viol=0;
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
             
          
             
            
             
             
             
             
             
             
             
             
             
             
             
             
         end
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         
         
         
%          mov{i_v}(k_mov(i_v)) = getframe(gcf);
%          k_mov(i_v)=k_mov(i_v)+1;
         
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
             
             x_fin=x_coord(served_nodes{i_l}(end-1)+1);
             y_fin=y_coord(served_nodes{i_l}(end-1)+1);
             served_node_coord=[x_fin,y_fin];
             dist_final=sqrt((served_node_coord(1)-x_coord(served_nodes{i_l}(end-2),1))^2+(served_node_coord(2)-y_coord(served_nodes{i_l}(end-2),1))^2)/v_speed(i_l,1);
             load_vec{i_l}=[load_vec{i_l};load_vec{i_l}(end)+q_vec(served_nodes{i_l}(end-1)+1,1)];
             start_time_mat{i_l}(end+1)=start_time_mat{i_l}(end)+dist_final+d_vec(served_nodes{i_l}(end-1)+1);
             
             current_time=start_time_mat{i_l}(end);
             %              figure(i_l)
             %              title(['Vehicle ',num2str(i_l),', Current time= ',num2str(current_time),', Global time= ',num2str(global_time)]);
             %              x_line=[x_coord(served_nodes{i_l}(end-1)+1) ,x_coord(served_nodes{i_l}(end)+1) ]';
             %              y_line=[y_coord(served_nodes{i_l}(end-1)+1) ,y_coord(served_nodes{i_l}(end)+1) ]';
             %
             %
             %              line(x_line,y_line,'Color','c','LineWidth',1);
             %              mov{i_l}(end) = getframe(gcf);
             %              filename=['anim' num2str(i_l) '.avi'];
             %              cd ('D:\VRP\Matlabstuff\codes\myAlgo-Dynamic\anim_multi')
             %              movie2avi(mov{i_l},filename,'fps',3, 'compression', 'None');
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
         
%         figure(i_l)
%          title(['Vehicle ',num2str(i_l),', Current time= ',num2str(current_time),', Global time= ',num2str(global_time)]);
%          x_line=[x_coord(served_nodes{i_l}(end-1)+1) ,x_coord(served_nodes{i_l}(end)+1) ]';
%          y_line=[y_coord(served_nodes{i_l}(end-1)+1) ,y_coord(served_nodes{i_l}(end)+1) ]';
         
         
        % line(x_line,y_line,'Color','m','LineWidth',1);
         %%%%animation creation%%%%%
%          mov{i_l}(end) = getframe(gcf);
%          filename=['anim' num2str(i_l) '.avi'];
%          cd ('/media/kaveh/DATA/VRP/Matlabstuff/codes/myAlgo-Dynamic (linux)/anim_multi')
%          movie2avi(mov{i_l},filename,'fps',3, 'compression', 'None');
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     else
         while served_nodes{i_l}(end)==served_nodes{i_l}(end-1)
             served_nodes{i_l}(end)=[];
             start_time_mat{i_l}(end)=[];
             load_vec{i_l}(end)=[];
         end
     end
     
     
     %%%%%%%% Adding routes to figures%%
     
     
     
     
     
     
     %%%%%Calculating the ride  time and excess ride time%%%%
     [ride_time{i_l},ride_viol(i_l)]=ride_time_cal(served_nodes{i_l},start_time_mat{i_l},d_vec,Max_ride);
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end
 
 
 x_d=[x_coordp(n_s) x_coordp(end-1)];
 y_d=[y_coordp(n_s) y_coordp(end-1)];

         point_coordinate_dyn=[x_d; y_d];
 for idy=1:m
       figure(idy)
        hold on
     scatter(point_coordinate_dyn(1,1),point_coordinate_dyn(2,1),'s','MarkerEdgeColor','r','MarkerFaceColor','blue');
             
              scatter(point_coordinate_dyn(1,2),point_coordinate_dyn(2,2),'s','MarkerEdgeColor','r','MarkerFaceColor','c');
              
              text(point_coordinate_dyn(1,1)+.5*rand,point_coordinate_dyn(2,1)+.4*rand,'d');
              
              text(point_coordinate_dyn(1,2)+.5*rand,point_coordinate_dyn(2,2)+.4*rand,'dd');
     
 
    
 end

 for i=1:m
%      x_coord=[];
%      y_coord=[];
     
      
      
      x_coord=[initial_depots_spec(i,1);x_coordp];
      
         y_coord=[initial_depots_spec(i,2);y_coordp];
     point_coordinate=[x_coord y_coord]';
    
     figure(i)
      hold on
           
              
     scatter(point_coordinate(1,2:n),point_coordinate(2,2:n),'o','MarkerEdgeColor','r','MarkerFaceColor','g');
     
     scatter(point_coordinate(1,n+2:2*n),point_coordinate(2,n+2:2*n),'o','MarkerEdgeColor','r','MarkerFaceColor','y');
     scatter(point_coordinate(1,1),point_coordinate(2,1),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
     scatter(point_coordinate(1,2*n+2),point_coordinate(2,2*n+2),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
     for ip=1:2*n+2
         %%%%to avoid putting initial depot no. and serving node
         %%%%together%%%
%          if x_coord(1)==x_coord(2)
%          
%          else
         text(point_coordinate(1,ip)+.1*rand,point_coordinate(2,ip)+.1*rand,num2str(ip-1));
%          end
     end
 
    
     title(['Vehicle ',num2str(i)])
    
     
 end
     
       
     
     
           for i_v=1:m
               figure(i_v)
               for i_line=1:size(served_nodes{i_v},1)-1
                   hold on
                   
                   c_t=start_time_mat{i_v}(i_line);
                 
                   title(['Vehicle ',num2str(i_v),', Current time= ',num2str(c_t),', Global time= ',num2str(global_time)]);
                  
            
             
             
              %%%%%%To start from each initial depot we need to change the first ponit's coordinate%%%%
              
              if i_line>1
                       x_line=[x_coord(served_nodes{i_v}(i_line)+1) ,x_coord(served_nodes{i_v}(i_line+1)+1)]';
                       y_line=[y_coord(served_nodes{i_v}(i_line)+1) ,y_coord(served_nodes{i_v}(i_line+1)+1)]';
                       
              else
                       x_line=[initial_depots_spec(i_v,1) ,x_coord(served_nodes{i_v}(i_line+1)+1)]';
                       y_line=[initial_depots_spec(i_v,2) ,y_coord(served_nodes{i_v}(i_line+1)+1)]';
              end       
                  line(x_line,y_line,'Color','m','LineWidth',1);
               end
               
               if size(infeas_nodes_mat{i_v},1)>0
                   for i_in=1:size(infeas_nodes_mat{i_v},1)
                       x_infes=x_coord( infeas_nodes_mat{i_v}(i_in)+1);
                       y_infes=y_coord(infeas_nodes_mat{i_v}(i_in)+1);
                       scatter(x_infes,y_infes,'o','MarkerEdgeColor','r','MarkerFaceColor','r');
                   end
               end
               
               
           end
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
 
 
 
 
 tot_dist=zeros(m,1);
 for i=1:m
     tot_dist(i)=distance_travelled(x_coord,y_coord,served_nodes{i});
 end
 Total_distance_traveled=sum(tot_dist);
 beta=2*n/(size(served_node_all,1)-1);
 Total_dist_modified=-(1+beta^gamma*(beta-1))* Total_distance_traveled;
 %%%%%%% Loop for solving darp ends here%%%%%%%%%%%