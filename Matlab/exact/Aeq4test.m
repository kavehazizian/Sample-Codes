clc;
clearvars;
close all
n=2;
m=2;
arc_set=arcs_cal(n);
t_x=4*n^2-n;
Num_var=4*m*(n*n+n+1);
A_eq4=zeros(m*2*n,Num_var);
for k=1:m
    for i=1:2*n
        
        for s=1:t_x
        if   arc_set(s,1)==i
                 
            A_eq4(i+(k-1)*2*n,s+(k-1)*t_x)=1;
        elseif arc_set(s,2)==i
            A_eq4(i+(k-1)*2*n,s+(k-1)*t_x)=-1;
        end
        end
    end
end
