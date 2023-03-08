function feasible=feasibility_chek(current_time,x_cur,y_cur,x_1,x_2,y_1,y_2,l_1,l_2,speed)
d_1=sqrt((x_cur-x_1)^2+(y_cur-y_1)^2)/speed;
d_2=sqrt((x_1-x_2)^2+(y_1-y_2)^2)/speed;
if (current_time+d_1>l_1)||(current_time+d_1+d_2>l_2)
    
    feasible=0;
   
else
    feasible=1;
end