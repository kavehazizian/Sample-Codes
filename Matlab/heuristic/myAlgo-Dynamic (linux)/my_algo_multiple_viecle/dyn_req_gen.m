function [x_d, y_d, d_d, q_d,e_d,l_d]=dyn_req_gen(g_time,x_range,y_range,cap,max_route)
 alp_dy=rand;
 bet_dy=rand;
 alp_dy_dest=rand;
 bet_dy_dest=rand;
 alpha_t=rand;
 
 x_dyn=alp_dy*x_range(1)+(1-alp_dy)*x_range(2);
 y_dyn=bet_dy*y_range(1)+(1-bet_dy)*y_range(2);
  x_dyn_dest=alp_dy_dest*x_range(1)+(1-alp_dy_dest)*x_range(2);
 y_dyn_dest=bet_dy_dest*y_range(1)+(1-bet_dy_dest)*y_range(2);
 %tim_dyn=(1-alpha_t)*g_time+alpha_t*max_route;
 %%%% Quick demand for a vehicle
  tim_dyn=g_time;
 %%%%
 e_or=tim_dyn;
 l_or=e_or+60;
 e_des=0;
 l_des=1400;
 
 
 x_d=[x_dyn x_dyn_dest];
 y_d=[y_dyn y_dyn_dest];
 q_dyn=ceil(alp_dy*cap+(1-alp_dy));
 d_dyn=q_dyn;
 
 q_d=[q_dyn -q_dyn];
 d_d=[d_dyn d_dyn];
 e_d=[e_or e_des];
 l_d=[l_or l_des];
 