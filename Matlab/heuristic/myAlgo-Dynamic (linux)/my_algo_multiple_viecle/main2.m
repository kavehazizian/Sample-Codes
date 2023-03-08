
   
 
clc
close all
clear all
[nodes,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec] = textread('pr01.txt','%d %f %f %d %d %f %f');
[v_no,v_cap,v_routing_time,v_speed] = textread('v2.txt','%d %d %f %f');
 n_s=(size(q_vec,1)-2)/2;
 m=size(v_no,1);
 max_route=max(v_routing_time);
  Max_ride=90;
 %%%%%% For dynamic part%%%%%
 load bestresult2;
 req_time=230;
 [initial_depots_spec,served_nodesd,start_time_matd,infeas_nodes_matd,load_vecd,current_time_matd,waiting_time_matd,...
x_coordd,y_coordd,d_vecd, q_vecd,e_vecd, l_vecd]=dyn_mod(req_time,served_nodes,start_time_mat,waiting_time_mat,...
load_vec,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec);
 
 %%%%%%%%%%%%%%%%%%%%%%%
 
 

%%%%%PSO Algo%%%%%
 pop_size=10;
 i_max=40;
 c_1=1.5;
 c_2=2;
 
 alpha=randi(20,pop_size,6);


 vel_vec=zeros(pop_size,6);
 p_best_vec=alpha;
 g_best_val=-inf;
 tet_min=.4;
 tet_max=4.9;
 init_time=req_time*ones(m,pop_size);
 
%%%%%%%%%%%%%%%%%%
i_pso=1;
i_chek=2;

while i_pso<=i_max && (i_chek-1)~=pop_size
    r_1=rand;
    r_2=rand;
     teta=tet_max-((tet_max-tet_min)/i_max)*i_pso;
    
     for i_pop=1:pop_size
         
         
         %%%%%%% Loop for solving darp starts here%%%%%%%%%%%
         [next_init(1:m,i_pop),Total_dist_modified_mat(i_pso,i_pop)]=core_part2(alpha(i_pop,:),x_coordd,y_coordd,d_vecd,...
             q_vecd,e_vecd, l_vecd,m,max_route, Max_ride,v_cap,v_routing_time,...
             v_speed,req_time,initial_depots_spec,served_nodesd,start_time_matd,load_vecd,waiting_time_matd);
         
         if i_pso>1&&(Total_dist_modified_mat(i_pso,i_pop)>Total_dist_modified_mat(i_pso-1,i_pop))
             p_best_vec(i_pop,:)=alpha(i_pop,:);
             
         end
         
         
         if Total_dist_modified_mat(i_pso,i_pop)>g_best_val
             
             g_best=alpha(i_pop,:);
             g_best_val=Total_dist_modified_mat(i_pso,i_pop);
             
             init_time_last=init_time;
         end
         
           
     end
     
     for i_pop2=1:pop_size
         
         vel_vec(i_pop2,:)=teta*vel_vec(i_pop2,:)+c_1*r_1*(p_best_vec(i_pop2,:)-alpha(i_pop2,:))+c_2*r_2*(g_best-alpha(i_pop2,:));
         
         %%%%%%%%Keeping the alphas inside the boundries, 0<alpha<100%%%%
         for i_bound=1:6
            
             if vel_vec(i_pop2,i_bound)>0
             vel_vec(i_pop2,i_bound)=vel_vec(i_pop2,i_bound)/(1+abs(vel_vec(i_pop2,i_bound)/(100-alpha(i_pop2,i_bound))));
             else
             vel_vec(i_pop2,i_bound)=vel_vec(i_pop2,i_bound)/(1+abs(vel_vec(i_pop2,i_bound)/(alpha(i_pop2,i_bound)-0)));   
                 
             end
             
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         alpha(i_pop2,:)=alpha(i_pop2,:)+vel_vec(i_pop2,:);
     end
    init_time=next_init; 
    %%%%%%%Checking the convergence%%%%%%
delta_T=0;    
i_chek=2;
epsil=1;
while delta_T<epsil&&i_chek<=pop_size
    %%%%%Stopping based on position of alpha%%%
    %delta_T=sum(abs(alpha(i_chek,:)-alpha(i_chek-1,:)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%Stopping based on objective function value%%%%
    delta_T=abs(Total_dist_modified_mat(i_pso,i_chek)-Total_dist_modified_mat(i_pso,i_chek-1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i_chek=i_chek+1;
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    i_pso=i_pso+1;
    
end

 
[ Total_distance_traveled,Total_dist_modified,served_nodesfinal,start_time_mat_final,infeas_nodes_mat_final,load_vec_final,...
    waiting_time_mat_final, ride_time_final]=final_result2(g_best,...
    x_coordd,y_coordd,d_vecd, q_vecd,e_vecd,l_vecd,m,max_route, Max_ride,v_cap,v_routing_time,v_speed,req_time,...
    initial_depots_spec,served_nodesd,start_time_matd,load_vecd,waiting_time_matd);