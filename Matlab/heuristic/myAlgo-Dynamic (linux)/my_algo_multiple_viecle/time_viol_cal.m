function t_v=time_viol_cal(served_nodes,start_time_mat,l_vec)
t_v_vec=zeros(size(served_nodes,1),1);
size(served_nodes,1);


for i=1:size(served_nodes,1)
    id= served_nodes(i,1);


 
    difference_time=start_time_mat(i)-l_vec(id+1,1);
    if difference_time>0
        t_v_vec(i,1)=difference_time;
    end
end
t_v=sum(t_v_vec(:,1));