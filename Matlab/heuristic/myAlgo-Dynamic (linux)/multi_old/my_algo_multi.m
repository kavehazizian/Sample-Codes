clc
close all
clear all
[nodes,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec] = textread('b2_16_small.dat','%d %f %f %d %d %f %f');
 n=(size(q_vec,1)-2)/2;
 %%%%%Modify this part for different Examples
 m=1;
 Max_ride=45;
 Cap_vec=6*ones(m,1);
  T_vec=480*ones(m,1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  

  
  
  current_time=0;
  
  cap_leaving=0;
  
  
  
  serving_node_test=0;
  served_nodes=0;
  i_nearet=1;
  

 
  
  
  
  %%%%%%%%%%%%%%%%%%%%Static Nodes%%%%%%%%%%%%%%%%%%%%%
  for i_f=1:m
      point_coordinate=[x_coord y_coord]';
      figure(i_f)
      scatter(point_coordinate(1,2:n+1),point_coordinate(2,2:n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','g');
      hold on
      scatter(point_coordinate(1,n+2:2*n+1),point_coordinate(2,n+2:2*n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','y');
      scatter(point_coordinate(1,1),point_coordinate(2,1),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
      scatter(point_coordinate(1,2*n+2),point_coordinate(2,2*n+2),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
      for ip=1:2*n+2
          text(point_coordinate(1,ip)+.1*rand,point_coordinate(2,ip)+.1*rand,num2str(ip-1));
      end
  end
  hold on
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  n_st=(size(e_vec,1)-2)/2;
  size_dyn=0;
  i_dyn=0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  no_neighbor=0;
  
  
  
  
  
  current_time_vec=zeros(1,m);
  cap_leaving_vec=zeros(1,m);
  speed_vec=ones(m,1);
  max_route_time=max(T_vec);
  
  served_node_coord_mat(1:2,1:m)=zeros(2,m);
  
  infeas_list=[];
  i_wait_vec=ones(m,1);
  waiting_time_mat=[];
  serving_node_test=0;
  served_nodes=0;
  remaining_v=ones(m,1);
  for i=1:m
      served_nodes_v{i}=0;
      Start_service{i}=0;
      load_leaving{i}=0;
      waiting_time_mat{i}=zeros(1,2);
      infeas_list_v{i}=[];
      i_time{i}=1;
  end
  while (min(current_time_vec)<=max_route_time&&sum(remaining_v)>0)
                  alp_dy=rand;
                  bet_dy=rand;
                  alp_dy_dest=rand;
                  bet_dy_dest=rand;
                  alpha_t=rand;
      for i_v=1:m
          v_cap=Cap_vec(i_v,1);
          current_time=current_time_vec(i_v);
          speed=speed_vec(i_v);
          cap_leaving=cap_leaving_vec(i_v);
          served_node_coord=served_node_coord_mat(1:2,i_v);
          
          i_wait= i_wait_vec(i_v,1);
          if current_time<T_vec(i_v)
              
              
              served_current_vh=served_nodes_v{i_v};
              infeas_list=infeas_list_v{i_v};
              [remaining,accepted_node,infeasible_node,update_time,leaving_load,waiting_time]=ddarp_one(x_coord,y_coord,e_vec,l_vec,q_vec,v_cap,Max_ride,served_node_coord,served_nodes,served_current_vh,current_time,speed,d_vec,infeas_list,cap_leaving);
              
              remaining_v(i_v,1)=remaining;
              infeas_list_v{i_v}=[infeas_list_v{i_v};infeasible_node];
              
              
              
              served_nodes=[served_nodes,accepted_node];
              i_time{i_v}=i_time{i_v}+1;
              
              
              served_nodes_v{i_v}=[served_nodes_v{i_v};accepted_node];
              waiting_time_mat{i_v}=[waiting_time_mat{i_v};waiting_time];
              cap_leaving_vec(i_v)=leaving_load;
              
              load_leaving{i_v}(i_time{i_v})=cap_leaving_vec(i_v);
              
              
              served_node_coord_mat(1:2,i_v)=[x_coord( accepted_node+1,1),y_coord(accepted_node+1,1)]';
              
              current_time_vec(i_v)=update_time;
              
              Start_service{i_v}(i_time{i_v})=current_time_vec(i_v);
              delay_anim=1;
              
              
              
              
              
              nod_plot_1=served_nodes_v{i_v}(end-1);
              nod_plot_2=served_nodes_v{i_v}(end);
              
              
              figure(i_v)
             
              title(['Current time= ',num2str(current_time)]);
              x_line=[x_coord(nod_plot_1+1) ,x_coord(nod_plot_2+1) ]';
              y_line=[y_coord(nod_plot_1+1) ,y_coord(nod_plot_2+1) ]';
              
              
              line(x_line,y_line,'Color','m','LineWidth',1);
              
              
              for i_nfe=1:size(infeas_list_v{i_v},1)
                  x_infes=x_coord(infeas_list_v{i_v}(i_nfe)+1);
                  y_infes=y_coord(infeas_list_v{i_v}(i_nfe)+1);
                  scatter(x_infes,y_infes,'o','MarkerEdgeColor','r','MarkerFaceColor','r');
              end
              pause(delay_anim)
              
              
              
              
              %%%%%%%Adding dynamic part and updating vectors%%%%%%%%%%%
              if i_dyn<size_dyn
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   alp_dy=rand;
%                   bet_dy=rand;
%                   alp_dy_dest=rand;
%                   bet_dy_dest=rand;
%                   alpha_t=rand;
                  tim_dyn= (1-alpha_t)*current_time+alpha_t*max_route_time;
                  x_dyn=alp_dy*(-10)+(1-alp_dy)*10;
                  y_dyn=bet_dy*(-10)+(1-bet_dy)*10;
                  x_dyn_dest=alp_dy_dest*(-10)+(1-alp_dy_dest)*10;
                  y_dyn_dest=bet_dy_dest*(-10)+(1-bet_dy_dest)*10;
                  e_or=tim_dyn;
                  l_or=e_or+rand*(max_route_time-e_or);
                  %  e_dest=e_or+Max_ride+rand*(max_route_time-e_or)
                  %  l_dest=e_dest+Max_ride+rand*abs((max_route_time-l_or))
                  e_dest=0;
                  l_dest=1440;
                  q_dyn=ceil(alp_dy*v_cap+(1-alp_dy));
                  d_dyn=q_dyn;
                  point_coordinate_dyn=[x_dyn x_dyn_dest;y_dyn y_dyn_dest];
                  
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                  figure(i_v)
                  
                  scatter(point_coordinate_dyn(1,1),point_coordinate_dyn(2,1),'s','MarkerEdgeColor','r','MarkerFaceColor','blue');
                  scatter(point_coordinate_dyn(1,2),point_coordinate_dyn(2,2),'s','MarkerEdgeColor','r','MarkerFaceColor','c');
                  
                  text(point_coordinate_dyn(1,1)+.1*rand,point_coordinate_dyn(2,1)+.1*rand,{'d_',num2str(i_dyn+1)});
                  text(point_coordinate_dyn(1,2)+.1*rand,point_coordinate_dyn(2,2)+.1*rand,{'dd_',num2str(i_dyn+1)});
                  
                  
                  x_coord=[x_coord(1:n+1);x_dyn;x_coord(n+2:2*n+1);x_dyn_dest;x_coord(end)];
                  y_coord=[y_coord(1:n+1);y_dyn;y_coord(n+2:2*n+1);y_dyn_dest;y_coord(end)];
                  d_vec=[d_vec(1:n+1);d_dyn;d_vec(n+2:2*n+1);d_dyn;d_vec(end)];
                  q_vec=[q_vec(1:n+1);q_dyn;q_vec(n+2:2*n+1);-q_dyn;q_vec(end)];
                  e_vec=[e_vec(1:n+1);e_or;e_vec(n+2:2*n+1);e_dest;e_vec(end)];
                  l_vec=[l_vec(1:n+1);l_or;l_vec(n+2:2*n+1);l_dest;l_vec(end)];
                  %  %%%%%%%%%%%%%%%%%%%%
                  n=(size(e_vec,1)-2)/2;
                 
                  
                  %%%%%%%%%%%%%Updating the number of seved nodes and infeasible nodes%%%
                 
                  for i_reord=1:size(served_nodes_v{i_v},1)
                      if served_nodes_v{i_v}(i_reord,1)>=n
                          served_nodes_v{i_v}(i_reord,1)=served_nodes(i_reord,1)+1 ;
                      end
                  end
                  
                  for i_reordp=1:size(infeas_list_v{i_v},1)
                      if infeas_list_v{i_v}(i_reordp,1)>=n
                          infeas_list_v{i_v}(i_reordp,1)=infeas_list_v{i_v}(i_reordp,1)+1;
                      end
                  end
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                  
                  
              end
              
          end
          
      end
      i_dyn=i_dyn+1; 
  end
  
  
  
  
 %%%%%%%%%%%%%%%%%%%last line to deopt%%%%%%%%%%%%%%%%%%%
 for i_l=1:m
 if served_nodes_v{i_l}(end)~=2*n+1;
served_nodes_v{i_l}=[served_nodes_v{i_l};2*n+1];
x_fin=x_coord(served_nodes_v{i_l}(end-1)+1);
y_fin=y_coord(served_nodes_v{i_l}(end-1)+1);
served_node_coord=[x_fin,y_fin];
dist_final=sqrt((served_node_coord(1)-x_coord(2*n+2,1))^2+(served_node_coord(2)-y_coord(2*n+2,1))^2)/speed_vec(i_l,1);
load_leaving{i_l}=[load_leaving{i_l},0];
Start_service{i_l}(end+1)=Start_service{i_l}(end)+dist_final+d_vec(served_nodes_v{i_l}(end-1)+1);

current_time=Start_service{i_l}(end);
     figure(i_l)
     title(['Current time= ',num2str(current_time)]);
     x_line=[x_coord(served_nodes_v{i_l}(end-1)+1) ,x_coord(served_nodes_v{i_l}(end)+1) ]';
     y_line=[y_coord(served_nodes_v{i_l}(end-1)+1) ,y_coord(served_nodes_v{i_l}(end)+1) ]';
     
     
     line(x_line,y_line,'Color','m','LineWidth',1);
 else
     while served_nodes_v{i_l}(end)==served_nodes_v{i_l}(end-1)
     served_nodes_v{i_l}(end)=[];
     Start_service{i_l}(end)=[];
     load_leaving{i_l}(end)=[];
     end
 end      
 end
 
 
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
 
 
 
 