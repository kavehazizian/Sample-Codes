function up_no=updating_no(x_vec,y_vec,d_vec, x_vecd,y_vecdd,d_vecd,node_no)

i=1;

while x_vec(node_no+1)~= x_vecd(i)||y_vec(node_no+1)~= y_vecdd(i)  || d_vec(node_no+1)~= d_vecd(i)
  
    i=i+1 ;
    
end

up_no=i;