function tot_dist=distance_travelled(x_coord,y_coord,served_nodes)

tot_dist=0;
 s_s=size(served_nodes,1);
 i=1;
 while i<s_s
     x_1=x_coord(served_nodes(i,1)+1,1);
     x_2=x_coord(served_nodes(i+1,1)+1,1);
     y_1=y_coord(served_nodes(i,1)+1,1);
     y_2=y_coord(served_nodes(i+1,1)+1,1);
    tot_dist=tot_dist+sqrt((x_1-x_2)^2+(y_1-y_2)^2); 
    i=i+1;
 end