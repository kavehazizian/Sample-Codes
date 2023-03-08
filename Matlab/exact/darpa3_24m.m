 
clc
clearvars
close all
[nodes,x_coord,y_coord,d_vec, q_vec,e_vec, l_vec] = textread('darpa3_24.dat','%d %f %f %d %d %f %f');
 m=3;
 n=30;
 Max_ride=45;
 Cap_vec=6*ones(m,1);
  T_vec=480*ones(m,1);
  arc_set=arcs_cal(n);
  t_x=4*n^2-n;
 t_vec_vec=zeros(t_x,1);
  %t_vec=zeros(t_x,1);
  for i=1:t_x
     ip=arc_set(i,1); 
     jp=arc_set(i,2);
      t_vec(i,1)=sqrt((x_coord(ip+1,1)-x_coord(jp+1,1))^2+(y_coord(ip+1,1)-y_coord(jp+1,1))^2);
    
  end
cost_matrix_vec=zeros(t_x*m,1);
for k=1:m
    cost_matrix_vec(1+(k-1)*t_x:k*t_x,1)=t_vec;
    
    
end