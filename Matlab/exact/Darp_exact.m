clc
clearvars
close all
m=2;
n=5;

Num_var=4*m*(n*n+n+1);

%%%%%Objective function
%input_data;
cost_matrix_vec=40*ones((4*n*n-n)*m,1);
f_vec=zeros(Num_var,1);

f_vec(1:(4*n*n-n)*m,1)=cost_matrix_vec;

%%%%%% Equality constraints
A_eq1=zeros(n,Num_var);
b_eq1=ones(n,1);
t_x=4*n*n-n;%%%vector x dimension
for i=1:n
    for j=1:m
    A_eq1(i,(2*n-1)*(i-1)+n+1+(j-1)*t_x:n+i*(2*n-1)+(j-1)*t_x)=ones(1,2*n-1);
    end
end

A_eq2=zeros(2*m,Num_var);
b_eq2=ones(2*m,1);
t_n=4*n*n-2*n;%%% starting elements for delivery arcs
for i=1:m
A_eq2(i,t_x*(i-1)+1:t_x*(i-1)+n)=ones(1,n);
A_eq2(m+i,t_x*(i-1)+t_n+1:t_x*(i-1)+t_n+n)=ones(1,n);
end

A_eq3=zeros(n*m,Num_var);
b_eq3=zeros(n*m,1);
for j=1:m
    for i=1:n
   A_eq3(i+(j-1)*n,(2*n-1)*(i-1)+n+1+(j-1)*t_x:n+i*(2*n-1)+(j-1)*t_x)=ones(1,2*n-1); 
   A_eq3(i+(j-1)*n,(2*n-1)*(i-1)+n+1+(j-1)*t_x+(2*n-1)*(n+1-i)+(i-1)*(2*n-2):n+i*(2*n-1)+(j-1)*t_x+(2*n-1)*(n+1-i)+(i-1)*(2*n-2)-1)=-ones(1,2*n-2);
   A_eq3(i+(j-1)*n,j*t_x-n+i)=-1; 
    end
end

A_eq4=zeros(2*n*m,Num_var);
b_eq4=zeros(2*n*m,1);

arc_set=arcs_cal(n);
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
% for j=1:m
%     for i=1:2*n
%         
%         
%         
%         if i<=n
%             alpha=n+i*(2*n-1)+(j-1)*t_x;
%             A_eq4(i+(j-1)*2*n,(2*n-1)*(i-1)+n+1+(j-1)*t_x:alpha)=ones(1,2*n-1);
%             % A_eq4(i+(j-1)*n,(j-1)*t_x+i)=-1;
%             
%             for k=1:2*n
%                 if k==1
%                     A_eq4(i+(j-1)*2*n,i+(j-1)*t_x)=-1;
%                 elseif k<=n
%                     if k<=i
%                         A_eq4(i+(j-1)*2*n,i+n+(2*n-1)*(i-k)-1+(j-1)*t_x)=-1;
%                     else
%                         A_eq4(i+(j-1)*2*n,i+n+(2*n-1)*(k-1)+(j-1)*t_x)=-1;
%                     end
%                 elseif k<n+i
%                     A_eq4(i+(j-1)*2*n,i+n+(2*n-1)*n+(k-n-1)*(2*n-2)-1+(j-1)*t_x)=-1;
%                 elseif k>n+i
%                     A_eq4(i+(j-1)*2*n,i+n+(2*n-1)*n+(k-n-1)*(2*n-2)+(j-1)*t_x)=-1;
%                 end
%                 
%                 
%             end
%         else
%             beta=n+n*(2*n-1)+(i-n)*(2*n-2)+(j-1)*t_x;
%             A_eq4(i+(j-1)*2*n,beta-(2*n-3):beta)=ones(1,2*n-2);
%             A_eq4(i+(j-1)*2*n,n+n*(2*n-1)+n*(2*n-2)+(i-n)+(j-1)*t_x)=1;
%             
%             for kp=1:2*n
%                 if kp<=n
%                     
%                     A_eq4(i+(j-1)*2*n,i+n+(2*n-1)*(kp-1)-1+(j-1)*t_x)=-1;
%                     
%                 elseif kp<i
%                     
%                     A_eq4(i+(j-1)*2*n,i+n+(2*n-1)*n+(kp-n-1)*(2*n-2)-2+(j-1)*t_x)=-1;
%                 elseif kp>i
%                     A_eq4(i+(j-1)*n,i+n+(2*n-1)*n+(kp-n-1)*(2*n-2)-1+(j-1)*t_x)=-1;
%                 end
%             end
%         end
%         
%         
%     end
% end

A_ineq1=zeros((4*n*n-n)*m, Num_var);
b_ineq1=zeros((4*n*n-n)*m,1);
A_ineq2=zeros((4*n*n-n)*m, Num_var);
b_ineq2=zeros((4*n*n-n)*m,1);

M_vec=1000*ones((4*n*n-n)*m,1);
W_vec=200*ones((4*n*n-n)*m,1);
t_vec=30*ones(4*n*n-n,1);
d_vec=3*ones(2*n+2,1);
d_vec(1,1)=0;
d_vec(2*n+2,1)=0;

q_vec=zeros(2*n+2,1);

q_vec(2:n+1,1)=ones(n,1);
q_vec(n+2:2*n+1,1)=-q_vec(2:n+1,1);

 counter=1;
 for k=1:m
     for j=1:2*n
         
         
         
         
         
         for i=0:2*n
             if j<=n
                 if i==0
                     A_ineq1(counter,j+(k-1)*(4*n*n-n))=M_vec(j+(k-1)*(4*n*n-n),1);
                     A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+1)=1;
                     A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+j+1)=-1;
                     b_ineq1(counter,1)=M_vec(j+(k-1)*(4*n*n-n),1)-d_vec(i+1,1)-t_vec(j,1);
                     
                     A_ineq2(counter,j+(k-1)*(4*n*n-n))=W_vec(j+(k-1)*(4*n*n-n),1);
                     A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+1)=1;
                     A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+j+1)=-1;
                     b_ineq2(counter,1)=W_vec(j+(k-1)*(4*n*n-n),1)-q_vec(j+1,1);
                     
                     
                     counter=counter+1;
                 elseif i<=n &&(j~=i&&i~=n+j)
                     cp=j;
                     for is=1:j
                     if i==is || i==is+n
                         cp=cp-1;
                     end
                     end
                     A_ineq1(counter,n+(i-1)*(2*n-1)+cp+(k-1)*(4*n*n-n))=M_vec((i-1)*(2*n-1)+j+(k-1)*(4*n*n-n),1);
                     A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+j+1)=-1;
                     A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+i+1)=1;
                     b_ineq1(counter,1)=M_vec((i-1)*(2*n-1)+j+(k-1)*(4*n*n-n),1)-d_vec(i+1,1)-t_vec((i-1)*(2*n-1)+j,1);
                     
                     A_ineq2(counter,n+(i-1)*(2*n-1)+cp+(k-1)*(4*n*n-n))=W_vec((i-1)*(2*n-1)+j+(k-1)*(4*n*n-n),1);
                     A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+j+1)=-1;
                     A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+i+1)=1;
                     b_ineq2(counter,1)=W_vec((i-1)*(2*n-1)+j+(k-1)*(4*n*n-n),1)-q_vec(j+1,1);
                     
                     
                     
                     counter=counter+1;
                 elseif i>n &&(j~=i&&i~=n+j)
                      cp=j;
                     for is=1:j
                     if i==is || i==is+n
                         cp=cp-1;
                         
                     end
                     end
                     A_ineq1(counter,n+n*(2*n-1)+(i-n-1)*(2*n-2)+cp+(k-1)*(4*n*n-n))=M_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j+(k-1)*(4*n*n-n),1);
                     A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+j+1)=-1;
                     A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+i+1)=1;
                     b_ineq1(counter,1)=M_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j+(k-1)*(4*n*n-n),1)-d_vec(i+1,1)-t_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j,1);
                     
                     A_ineq2(counter,n+n*(2*n-1)+(i-n-1)*(2*n-2)+cp+(k-1)*(4*n*n-n))=W_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j+(k-1)*(4*n*n-n),1);
                     A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+j+1)=-1;
                     A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+i+1)=1;
                     b_ineq2(counter,1)=W_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j+(k-1)*(4*n*n-n),1)-q_vec(j+1,1); 
                     
                     
                     counter=counter+1;
                 end
             elseif i<=n &&i~=j &&i~=0 
                 
                 cp=j;
                     for is=1:j
                     if i==is || i==is+n
                         cp=cp-1;
                         
                     end
                     end
                     
                 A_ineq1(counter,n+(i-1)*(2*n-1)+cp+(k-1)*(4*n*n-n))=M_vec((i-1)*(2*n-1)+j-1+(k-1)*(4*n*n-n),1);
                 A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+j+1)=-1;
                 A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+i+1)=1;
                 b_ineq1(counter,1)=M_vec((i-1)*(2*n-1)+j-1+(k-1)*(4*n*n-n),1)-d_vec(i+1,1)-t_vec((i-1)*(2*n-1)+j-1,1);
                 
                 A_ineq2(counter,n+(i-1)*(2*n-1)+cp+(k-1)*(4*n*n-n))=W_vec((i-1)*(2*n-1)+j-1+(k-1)*(4*n*n-n),1);
                 A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+j+1)=-1;
                 A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+i+1)=1;
                 b_ineq2(counter,1)=W_vec((i-1)*(2*n-1)+j-1+(k-1)*(4*n*n-n),1)-q_vec(j+1,1);
                 
                 counter=counter+1;
                 
             
             elseif i>n &&i~=j 
                 cp=j;
                     for is=1:j
                     if i==is || i==is+n
                         cp=cp-1;
                         
                     end
                     end
                 A_ineq1(counter,n+n*(2*n-1)+(i-n-1)*(2*n-2)+cp+(k-1)*(4*n*n-n))=M_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j-1+(k-1)*(4*n*n-n),1);
                 A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+j+1)=-1;
                 A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+i+1)=1;
                 b_ineq1(counter,1)=M_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j-1+(k-1)*(4*n*n-n),1)-d_vec(i+1,1)-t_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j-1,1);
                 
                 A_ineq2(counter,n+n*(2*n-1)+(i-n-1)*(2*n-2)+cp+(k-1)*(4*n*n-n))=W_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j-1+(k-1)*(4*n*n-n),1);
                 A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+j+1)=-1;
                 A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+i+1)=1;
                 b_ineq2(counter,1)=W_vec(n*(2*n-1)+(i-n-1)*(2*n-2)+j-1+(k-1)*(4*n*n-n),1)-q_vec(j+1,1);
                 
                 
                 
                 counter=counter+1;
             end
             
             
             
             
             
                          
             
         end
         
     end
     
     
     
     
     
     
     j=2*n+1;
     
     for ip=1:n
         
         A_ineq1(counter,n+(2*n-1)*n+(2*n-2)*n+ip+(k-1)*(4*n*n-n))=M_vec((2*n-1)*n+(2*n-2)*n+ip+(k-1)*(4*n*n-n),1);
         A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+j+1)=-1;
         A_ineq1(counter,m*(4*n*n-n)+(k-1)*(2*n+2)+ip+n+1)=1;
         b_ineq1(counter,1)=M_vec((2*n-1)*n+(2*n-2)*n+ip+(k-1)*(4*n*n-n),1)-d_vec(ip+n+1,1)-t_vec((2*n-1)*n+(2*n-2)*n+ip,1);
         
         A_ineq2(counter,n+(2*n-1)*n+(2*n-2)*n+ip+(k-1)*(4*n*n-n))=W_vec((2*n-1)*n+(2*n-2)*n+ip+(k-1)*(4*n*n-n),1);
         A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+j+1)=-1;
         A_ineq2(counter,m*(4*n*n-n)+m*(2*n+2)+(k-1)*(2*n+2)+ip+n+1)=1;
         b_ineq2(counter,1)=W_vec((2*n-1)*n+(2*n-2)*n+ip+(k-1)*(4*n*n-n),1)-q_vec(j+1,1);
         
         
         counter=counter+1;
     end
     
     
     
     
 end
     
A_ineq3=zeros(m*n,Num_var);
b_ineq3=zeros(n*m,1);
Num_ineq3=m*(4*n^2+3*n+4);

cout_3=1;
for k=1:m
    
    for i=1:n

 A_ineq3(cout_3,Num_ineq3+i+(k-1)*n)=-1;
  A_ineq3(cout_3,m*t_x+i+n+1+(k-1)*(2*n+2))=1;
  A_ineq3(cout_3,m*t_x+i+1+(k-1)*(2*n+2))=-1;
  b_ineq3(cout_3,1)=d_vec(i+1,1);
 cout_3=cout_3+1;
    end
end

T_vec=3600*ones(m,1);

A_ineq4=zeros(m,Num_var);
b_ineq4=T_vec;

for k=1:m
    
 A_ineq4(k,m*t_x+(k-1)*(2*n+2)+1)=-1;   
 A_ineq4(k,m*t_x+k*(2*n+2))=1;   
  
end


%%%%%% Bounds on Variables%%%
lb=zeros(Num_var,1);
ub=inf*ones(Num_var,1);
ub(1:t_x*m,1)=ones(t_x*m,1);
e_vec=0*ones(2*n+2,1);
l_vec=520*ones(2*n+2,1);
%%%%%changing time windows
e_vec(n+2:2*n+2,1)=1*ones(n+1,1);
l_vec(1:n+1,1)=510*ones(n+1,1);
%%%%%%%%%%%%%%%%%%%%%%%%

Max_ride=45;
Cap_vec=3*ones(m,1);

for k=1:m
    lb(t_x*m+(k-1)*(2*n+2)+1:t_x*m+k*(2*n+2),1)=e_vec;
    ub(t_x*m+(k-1)*(2*n+2)+1:t_x*m+k*(2*n+2),1)=l_vec;
    for i=1:n
        lb(Num_var-m*n+i+(k-1)*n,1)=t_vec(n+(i-1)*(2*n-1)+(n-1)+i,1);
        
    end
    
    ub(Num_var-m*n+1:Num_var,1)=Max_ride*ones(m*n,1);
    
    for ip=1:2*n+2
        lb(Num_var-m*n-m*(2*n+2)+(k-1)*(2*n+2)+ip,1)=max(0,q_vec(ip));
        ub(Num_var-m*n-m*(2*n+2)+(k-1)*(2*n+2)+ip,1)=min(Cap_vec(k),Cap_vec(k)+q_vec(ip));
        
    end
    
    
    
end

A_ineq=[A_ineq1;A_ineq2;A_ineq3;A_ineq4];
b_ineq=[b_ineq1;b_ineq2;b_ineq3;b_ineq4];

A_eq=[A_eq1;A_eq2;A_eq3;A_eq4];
b_eq=[b_eq1;b_eq2;b_eq3;b_eq4];

int_var=1:m*t_x;
options = optimoptions('intlinprog','Display','off');
[x,fval,exitflag,output]  = intlinprog(f_vec,int_var,A_ineq,b_ineq,A_eq,b_eq,lb,ub,options)

%%%%%%%%Post Processing%%%
if size(x,1)>1
    Passed_arcs=x(1:m*t_x,1);
    Service_start=x(m*t_x+1:m*t_x+m*(2*n+2),1);
    Loads=x(m*t_x+m*(2*n+2)+1:m*t_x+2*m*(2*n+2),1);
    Ride_time=x(Num_var-n*m+1:Num_var,1);
    arcs_mat=zeros(t_x,m);
    service_start_mat= zeros(2*n+2,m);
    load_mat= zeros(2*n+2,m);
    ride_time_mat=zeros(n,m);
    for k=1:m
        arcs_mat(1:t_x,k)=Passed_arcs((k-1)*t_x+1:k*t_x,1)
        service_start_mat(1:2*n+2,k)=Service_start((k-1)*(2*n+2)+1:k*(2*n+2),1)
        load_mat(1:2*n+2,k)=Loads((k-1)*(2*n+2)+1:k*(2*n+2),1)
        ride_time_mat(1:n,k)=Ride_time(1+(k-1)*n:k*n,1)
    end
    %%%%%%%Visualization
    
    point_coordinate=10*rand(2,2*n+2);
    figure(1)
    scatter(point_coordinate(1,2:n+1),point_coordinate(2,2:n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','g');
    hold on
    scatter(point_coordinate(1,n+2:2*n+1),point_coordinate(2,n+2:2*n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','y');
    scatter(point_coordinate(1,1),point_coordinate(2,1),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
    scatter(point_coordinate(1,2*n+2),point_coordinate(2,2*n+2),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
    title('All possible arcs')
    hold on
    yp=zeros(2,n+1);
    zp=zeros(2,n+1);
    arc_count=0;
    
    
    for i=1:(2*n+2)
        text(point_coordinate(1,i)+.1*rand,point_coordinate(2,i)+.1*rand,num2str(i-1));
        
        for j=1:(2*n+2)
            
            if (i==j)||(i==j+n)||(i==1&&j>n+1)|| (j==2*n+2&&i<=n+1)||(j==1&&i>1)||(i<n+1&&j==2*n+2)||(i==2*n+2)
                
                
                %%%% No arcs
                
            else
                
                
                dist=sqrt((point_coordinate(1,i)-point_coordinate(1,j))^2+(point_coordinate(2,i)-point_coordinate(2,j))^2);
                zp=dist/2*rand(2,1);
                
                
                
                
                
                xpp=[point_coordinate(1,i);point_coordinate(1,j);zp(1)];
                ypp=[point_coordinate(2,i);point_coordinate(2,j);zp(2)];
                %%%%% Interpolation for degree of three
                %up=dist/2*rand(2,1);
                %             xpp=[point_coordinate(1,i);point_coordinate(1,j);zp(1);up(1)];
                %             ypp=[point_coordinate(2,i);point_coordinate(2,j);zp(2);up(2)];
                pol_camp= polyfit(xpp,ypp,2);
                
                
                xValues = linspace(min(xpp(1:2)),max(xpp(1:2)));
                y1 = polyval(pol_camp,xValues);
                plot(xValues,y1)
               
                arc_count=arc_count+1;
                if j==2*n+2
                countp=arc_count+(2*n+1-i)*(2*n-2)+(i-n-1-1);
                xval_mat{countp}=xValues;
               y_mat{countp}=y1;
               arc_count=arc_count-1;
                else     
               xval_mat{arc_count}=xValues;
               y_mat{arc_count}=y1;
                end
            end
        end
        
    end
    
    
    
    for k=1:m
        
        
        figure(k+1)
        
    scatter(point_coordinate(1,2:n+1),point_coordinate(2,2:n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','g');
    hold on
    scatter(point_coordinate(1,n+2:2*n+1),point_coordinate(2,n+2:2*n+1),'o','MarkerEdgeColor','r','MarkerFaceColor','y');
    scatter(point_coordinate(1,1),point_coordinate(2,1),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
    scatter(point_coordinate(1,2*n+2),point_coordinate(2,2*n+2),'o','MarkerEdgeColor','r','MarkerFaceColor','black');
                             for ip=1:2*n+2
                              text(point_coordinate(1,ip)+.1*rand,point_coordinate(2,ip)+.1*rand,num2str(ip-1));
                             end
                             title(['Arcs passed by vehicle ',num2str(k)])
                             hold on
        for i=1:(4*n^2-n)
        
   hold on
                             if (.99999<=arcs_mat(i,k)&&arcs_mat(i,k)<=1.111111)
                                plot(xval_mat{i},y_mat{i})
        
       
                               % hold on
                             end
        
        end  
    end
    

end