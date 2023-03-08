
   
 
clc
close all
clear all
[nodes,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec] = textread('pr01_small.txt','%d %f %f %d %d %f %f');
[v_no,v_cap,v_routing_time,v_speed] = textread('vehicles.txt','%d %d %f %f');
 n_s=(size(q_vec,1)-2)/2;
 m=size(v_no,1);
 max_route=max(v_routing_time);
 
 
 
 Max_ride=90;
allowed_delay=15;
 alpha=.6;

    
 global_time=0;
 pool_size=2*n_s+2;
 time_res=5;
    for i=1:m
    served_nodes{i}=[0];
    start_time_mat{i}=[0];
    infeas_nodes_mat{i}=[];
     load_vec{i}=[0];
     current_time_mat{i}=[0];
     waiting_time_mat{i}=zeros(1,2);
    end
    
 served_node_all=[0];
 
 %%%%%%Static Figures%%%

 n=n_s;
 for i=1:m
     point_coordinate=[x_coord y_coord]';
    figure(i)
    scatter(point_coordinate(1,2:n+1),point_coordinate(2,2:n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','g');
    hold on
    scatter(point_coordinate(1,n+2:2*n+1),point_coordinate(2,n+2:2*n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','y');
    scatter(point_coordinate(1,1),point_coordinate(2,1),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
    scatter(point_coordinate(1,2*n+2),point_coordinate(2,2*n+2),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
    for ip=1:2*n+2
        text(point_coordinate(1,ip)+.1*rand,point_coordinate(2,ip)+.1*rand,num2str(ip-1));
    end

 title(['Vehicle ',num2str(i)])

     
 end
 
 %%%%%%%%Initiating some variables%%%%
 k_mov=ones(m,1);
 n_d=7;
  
 i_dyn=1;
 x_range=[-10,10];
 y_range=[-10,10];                                                                                                                                                                      
 n=n_s+n_d;
 n_add=n_s;

 %%%%%%%%%%%%%%%%%%%%%%%%
 

 while (global_time<max_route&&pool_size>0)
     %alpha=alpha+ alpha_increment; 
     for i_v=1:m
         x_line=[];
         y_line=[];                             %%%% current_time_mat{i_v}(1) allows to start serving even after time zero%%%%   
         if current_time_mat{i_v}(end)<=global_time+current_time_mat{i_v}(1)
             cap=v_cap(i_v);
             rout_time=v_routing_time(i_v); 
             speed=v_speed(i_v);
          
             [served_node,start_time,infeas_nodes,current_time,load,waiting]=solve_single(served_node_all,served_nodes{i_v},infeas_nodes_mat{i_v},start_time_mat{i_v},current_time_mat{i_v}(end),load_vec{i_v}(end),cap,speed,Max_ride,x_coord,y_coord,d_vec,q_vec,e_vec,l_vec);
            
             
             
             served_nodes{i_v}=[served_nodes{i_v};served_node];
             start_time_mat{i_v}=[start_time_mat{i_v};start_time];
             infeas_nodes_mat{i_v}=[infeas_nodes_mat{i_v};infeas_nodes];
             current_time_mat{i_v}=[current_time_mat{i_v};current_time];
             load_vec{i_v}=[load_vec{i_v};load];
             
             waiting_time_mat{i_v}=[waiting_time_mat{i_v};waiting];
             if served_node>0
                 served_node_all=[served_node_all;served_node];
             end
             
          
             
             %%%%%%%% Adding routes to figures%%
             
             figure(i_v)
             
             hold on
             %pause(3)
             
             c_t=current_time_mat{i_v}(end);
             title(['Vehicle ',num2str(i_v),', Current time= ',num2str(c_t),', Global time= ',num2str(global_time)]);
             
             
             %               fig_vec(i_v)=fig_vec(i_v)+1;
             %               i_fig= fig_vec(i_v);
             if size(served_nodes{i_v},1)>1
             x_line=[x_coord(served_nodes{i_v}(end-1)+1) ,x_coord(served_nodes{i_v}(end)+1)]';
             y_line=[y_coord(served_nodes{i_v}(end-1)+1) ,y_coord(served_nodes{i_v}(end)+1)]';
             line(x_line,y_line,'Color','m','LineWidth',1);
             x_line=[];
             y_line=[];
             end
             if size(infeas_nodes_mat{i_v},1)>0
                 for i_in=1:size(infeas_nodes_mat{i_v},1)
                     x_infes=x_coord( infeas_nodes_mat{i_v}(i_in)+1);
                     y_infes=y_coord(infeas_nodes_mat{i_v}(i_in)+1);
                     scatter(x_infes,y_infes,'o','MarkerEdgeColor','r','MarkerFaceColor','r');
                 end
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
             
             
             
             
             
             
             
             
             
             
             
         end
         
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         
         
         
         mov{i_v}(k_mov(i_v)) = getframe(gcf);
         k_mov(i_v)=k_mov(i_v)+1;
         
     end
%       served_nodes{1}'
% %             n_add
% %             i_dyn
% %             dyn_req
% disp('here')
%              pause(2)
     
     %figure(i_v)
     %%%%%%If wanna show all steps of global time%%%
     %        pause(1)
     %          c_t=current_time_mat{i_v}(end);
     %           title(['Vehicle ',num2str(i_v),', Current time= ',num2str(c_t),', Global time= ',num2str(global_time)]);
     %%%%%%%%%%%%%%%%%%%%%%%%%%
     %            pause(1)
     %          x_line=[x_coord(served_nodes{i_v}(end-1)+1) ,x_coord(served_nodes{i_v}(end)+1)]';
     %               y_line=[y_coord(served_nodes{i_v}(end-1)+1) ,y_coord(served_nodes{i_v}(end)+1)]';
     %               line(x_line,y_line,'Color','m','LineWidth',1);
     % hold on
     
     %%%%%%%%%%%%%%%%%%%%%Adding dynamic part%%%%%%%
%       alpha_t=rand;       
%       request_time=(1-alpha_t)*global_time+alpha_t*max_route;
request_time=100;

      if i_dyn<=n_d&&global_time>=request_time
          
          [x_d, y_d, d_d, q_d,e_d,l_d]=dyn_req_gen(global_time,x_range,y_range,cap,max_route);
          
          
          
          
          
          for id=1:m
              point_coordinate_dyn=[x_d; y_d];
              figure(id)
              scatter(point_coordinate_dyn(1,1),point_coordinate_dyn(2,1),'s','MarkerEdgeColor','r','MarkerFaceColor','blue');
              scatter(point_coordinate_dyn(1,2),point_coordinate_dyn(2,2),'s','MarkerEdgeColor','r','MarkerFaceColor','c');
              
              text(point_coordinate_dyn(1,1)+.1*rand,point_coordinate_dyn(2,1)+.1*rand,{'d_',num2str(i_dyn)});
              text(point_coordinate_dyn(1,2)+.1*rand,point_coordinate_dyn(2,2)+.1*rand,{'dd_',num2str(i_dyn)});
          end
          
          
          x_coord=vec_aug(x_coord,x_d,n_add);
          y_coord=vec_aug(y_coord,y_d,n_add);
          d_vec=vec_aug(d_vec,d_d,n_add);
          q_vec=vec_aug(q_vec,q_d,n_add);
          e_vec=vec_aug(e_vec,e_d,n_add);
          l_vec=vec_aug(l_vec,l_d,n_add);
          
          n_add=n_s+i_dyn;
          
          
          
          
          
          
          
          %%%%%%%%%%%%%Updating the number of served nodes and infeasible nodes%%%
          for i_vv=1:m
              for i_reord=1:size(served_nodes{i_vv},1)
                  
                  if served_nodes{i_vv}(i_reord,1)>=n_add
                      served_nodes{i_vv}(i_reord,1)=served_nodes{i_vv}(i_reord,1)+1;
                      
                      served_nodes{i_vv}(i_reord,1)
                      
                  end
              end
              
              for i_reordp=1:size(infeas_nodes_mat{i_vv},1)
                  if infeas_nodes_mat{i_vv}(i_reordp,1)>=n_add
                      infeas_nodes_mat{i_vv}(i_reordp,1)=infeas_nodes_mat{i_vv}(i_reordp,1)+1;
                  end
              end
              
              for i_reordpp=1:size(waiting_time_mat{i_vv},1)
                  if waiting_time_mat{i_vv}(i_reordpp,1)>=n_add
                      waiting_time_mat{i_vv}(i_reordpp,1)=waiting_time_mat{i_vv}(i_reordpp,1)+1;
                  end
              end
              
              
          end
          
          for i_reordppp=2:size(served_node_all,1)
              if served_node_all(i_reordppp,1)>=n_add
                  
                  served_node_all(i_reordppp,1)=served_node_all(i_reordppp,1)+1;
                  
              end
          end
          
          i_dyn=i_dyn+1;
          
          
      end
      
      
      
      
      
      
      global_time=global_time+time_res;
      pool_size=2*n_s+1+2*n_d-size(served_node_all,1);
      
 end
 
 
 %%%%%%%%%%%%%%%%%%%last line to deopt%%%%%%%%%%%%%%%%%%%
 %n=n_s;
 for i_l=1:m
     
     %%%%%%%%%%%%%%%%%%%%Allowing the last picked-up customer to reach iits
     %%%%%%%%%%%%%%%%%%%%destination%%%%
     i_rem_co=size(served_nodes{i_l},1);
     still_onboard=[];
     for i_pass=i_rem_co:-1:ceil(i_rem_co/2)
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
             figure(i_l)
             title(['Vehicle ',num2str(i_l),', Current time= ',num2str(current_time),', Global time= ',num2str(global_time)]);
             x_line=[x_coord(served_nodes{i_l}(end-1)+1) ,x_coord(served_nodes{i_l}(end)+1) ]';
             y_line=[y_coord(served_nodes{i_l}(end-1)+1) ,y_coord(served_nodes{i_l}(end)+1) ]';
             
             
             line(x_line,y_line,'Color','c','LineWidth',1);
             mov{i_l}(end) = getframe(gcf);
             filename=['anim' num2str(i_l) '.avi'];
             cd ('D:\VRP\Matlabstuff\codes\myAlgo-Dynamic\anim_multi')
             movie2avi(mov{i_l},filename,'fps',3, 'compression', 'None');
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
         figure(i_l)
         title(['Vehicle ',num2str(i_l),', Current time= ',num2str(current_time),', Global time= ',num2str(global_time)]);
         x_line=[x_coord(served_nodes{i_l}(end-1)+1) ,x_coord(served_nodes{i_l}(end)+1) ]';
         y_line=[y_coord(served_nodes{i_l}(end-1)+1) ,y_coord(served_nodes{i_l}(end)+1) ]';
         
         
         line(x_line,y_line,'Color','m','LineWidth',1);
         mov{i_l}(end) = getframe(gcf);
         filename=['anim' num2str(i_l) '.avi'];
         cd ('D:\VRP\Matlabstuff\codes\myAlgo-Dynamic\anim_multi')
         movie2avi(mov{i_l},filename,'fps',3, 'compression', 'None');
     else
         while served_nodes{i_l}(end)==served_nodes{i_l}(end-1)
             served_nodes{i_l}(end)=[];
             start_time_mat{i_l}(end)=[];
             load_vec{i_l}(end)=[];
         end
     end
     %%%%Calculating the ride  time%%%%
     ride_time{i_l}=ride_time_cal(served_nodes{i_l},start_time_mat{i_l},n,d_vec);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 end
 
 
 tot_dist=zeros(m,1);
 for i=1:m
     tot_dist(i)=distance_travelled(x_coord,y_coord,served_nodes{i});
 end
 Total_distance_traveled=sum(tot_dist)
 %%%%%%%changing alpha%%%%
%  distal(i_alpha,1)=alpha;
%  distal(i_alpha,2)=size(served_node_all,1);
%  distal(i_alpha,3)=Total_distance_traveled;
%  alpha=alpha+1/20;


