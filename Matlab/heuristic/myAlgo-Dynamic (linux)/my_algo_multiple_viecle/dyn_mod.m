function [initial_depots_spec,served_nodesd,start_time_matd,infeas_nodes_matd,load_vecd,current_time_matd,waiting_time_matd,...
x_coordd,y_coordd,d_vecd, q_vecd,e_vecd, l_vecd]=dyn_mod(req_time,served_nodes,start_time_mat,waiting_time_mat,...
load_vec,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec)


m=size (served_nodes,2);
n_s=size(e_vec,1)/2-1;
rows2remove=[0];


initial_depots_spec=zeros(m,6);
for i=1:m
already_served_nodes{i}=[0];
incomplete_req{i}=[];
start_time_matd{i}=[];
waiting_time_matd{i}=[];
load_vecd{i}=[];
j=2;
while start_time_mat{i}(j,1)<req_time&&j<=size(served_nodes{i},1)
j=j+1;
end
already_served_nodes{i}=served_nodes{i}(1:j-1);
end
for i=1:m
for k=1:size(already_served_nodes{i},1)
if already_served_nodes{i}(k)>n_s
rows2remove=[rows2remove;already_served_nodes{i}(k)-n_s;already_served_nodes{i}(k)];
elseif ismember(already_served_nodes{i}(k)+n_s,already_served_nodes{i})==0&&already_served_nodes{i}(k)~=0
incomplete_req{i}=[incomplete_req{i};already_served_nodes{i}(k)];
end
end
initial_depots(i)=already_served_nodes{i}(end);
end

for ii=1:m
  
initial_depots_spec(ii,1:6)=[x_coord(initial_depots(ii)+1) y_coord(initial_depots(ii)+1) d_vec(initial_depots(ii)+1) q_vec(initial_depots(ii)+1) e_vec(initial_depots(ii)+1) l_vec(initial_depots(ii)+1)];

end
%%%%%%%%%%%%% Since row numbering in input files starts with 1 say node
%%%%%%%%%%%%% zer is in row one hence we must increase removing index by
%%%%%%%%%%%%% one
%%%%Keeping required static info for node no. updating purpose%%%
x_coordp=x_coord;
y_coordp=y_coord;
d_vecp=d_vec;
%%%%%%%%%%%%%%%
rows2remove=rows2remove+1;
%%%%deleting the already served nodes from  the candidate pool%%%
x_coord(rows2remove,:)=[];
y_coord(rows2remove,:)=[];
d_vec(rows2remove,:)=[];
q_vec(rows2remove,:)=[];
e_vec(rows2remove,:)=[];
l_vec(rows2remove,:)=[];
[x_d,y_d,d_d,q_d,e_d,l_d]=textread('dyn_input.txt','%f %f %d %d %f %f');

n_add=(size(e_vec,1)-1)/2;
x_coordd=vec_aug_dyn(x_coord,x_d,n_add);
y_coordd=vec_aug_dyn(y_coord,y_d,n_add);
d_vecd=vec_aug_dyn(d_vec,d_d,n_add);
q_vecd=vec_aug_dyn(q_vec,q_d,n_add);
e_vecd=vec_aug_dyn(e_vec,e_d,n_add);
l_vecd=vec_aug_dyn(l_vec,l_d,n_add);
for i=1:m
    %%%%Updating no. of incomplete served nodes in new pool%%%
    
    for i_m=1:size(incomplete_req{i},1)
        imm=1;
        %%%%%%Picking-up the start-time for  picked up
        %%%%%%customers%%
        [~,wr]=ismember(incomplete_req{i}(i_m),served_nodes{i});
        start_time_matd{i}(i_m,1)=start_time_mat{i}(wr);
        
        load_vecd{i}(i_m,1)=load_vec{i}(wr);
        [we,wrw]=ismember(incomplete_req{i}(i_m),waiting_time_mat{i}(:,1));
        if we==1
            waiting_time_matd{i}(imm,:)= waiting_time_mat{i}(wrw,:);
            imm=imm+1;
        end
        
        %%%%%%%%%
        
        up_no=updating_no(x_coordp,y_coordp,d_vecp, x_coordd,y_coordd,d_vecd,incomplete_req{i}(i_m));
        incomplete_req{i}(i_m,1)=up_no;
    end
    %%%Adding last served node as initial depot
    [~,wrp]=ismember(already_served_nodes{i}(end),served_nodes{i});
    
    start_time_matd{i}=[start_time_mat{i}(wrp);start_time_matd{i}];
    load_vecd{i}=[load_vec{i}(wrp);load_vecd{i}];
    served_nodesd{i}= [0;incomplete_req{i}];
    [wep,wrwp]=ismember(already_served_nodes{i}(end),waiting_time_mat{i}(:,1));
    if wep==1
        waiting_time_matd{i}=[waiting_time_mat{i}(wrwp,:);waiting_time_matd{i}] ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    j_break=size(already_served_nodes{i},1);
    
    infeas_nodes_matd{i}=[];
    
    current_time_matd{i}=req_time;
    
    
end

