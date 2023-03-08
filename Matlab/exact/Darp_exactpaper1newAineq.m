%inputdarp1
%pr01_small
inputdarp;
%inputdarp_small;
%darpb_16;
%darpa2_20;
%darpa4_40;
 %darpa3_24;
 %darpa3_36;
 %darpb3_36;
%darpa3_30;
%darpa8_96;
%%%%%%%Inputs%%%%%%
% m=2;
% n=3;
% cost_matrix_vec=40*ones((4*n*n-n)*m,1);
% M_vec=12000*ones((4*n*n-n)*m,1);
% W_vec=200*ones((4*n*n-n)*m,1);
% t_vec=30*ones(4*n*n-n,1);
% d_vec=3*ones(2*n+2,1);
% d_vec(1,1)=0;
% d_vec(2*n+2,1)=0;
% e_vec=0*ones(2*n+2,1);
% l_vec=520*ones(2*n+2,1);
% %%%%%changing time windows
% e_vec(n+2:2*n+2,1)=1*ones(n+1,1);
% l_vec(1:n+1,1)=510*ones(n+1,1);
% %%%%%%%%%%%%%%%%%%%%%%%%
% Max_ride=45;
% Cap_vec=3*ones(m,1);
% q_vec=sparse(2*n+2,1);
% 
% q_vec(2:n+1,1)=ones(n,1);
% q_vec(n+2:2*n+1,1)=-q_vec(2:n+1,1);
% 
% T_vec=3600*ones(m,1);


%%%%%%%%%%%%%%%%%%%
Num_var=4*m*(n*n+n+1);

%%%%%Objective function

f_vec=zeros(Num_var,1);
f_vec(1:(4*n*n-n)*m,1)=cost_matrix_vec;

%%%%%% Equality constraints
A_eq1=sparse(n,Num_var);
b_eq1=ones(n,1);
t_x=4*n*n-n;%%%vector x dimension
for i=1:n
    for j=1:m
    A_eq1(i,(2*n-1)*(i-1)+n+1+(j-1)*t_x:n+i*(2*n-1)+(j-1)*t_x)=ones(1,2*n-1);
    end
end

A_eq2=sparse(2*m,Num_var);
b_eq2=ones(2*m,1);
t_n=4*n*n-2*n;%%% starting elements for delivery arcs
for i=1:m
A_eq2(i,t_x*(i-1)+1:t_x*(i-1)+n)=ones(1,n);
A_eq2(m+i,t_x*(i-1)+t_n+1:t_x*(i-1)+t_n+n)=ones(1,n);
end

A_eq3=sparse(n*m,Num_var);
b_eq3=sparse(n*m,1);
for j=1:m
    for i=1:n
   A_eq3(i+(j-1)*n,(2*n-1)*(i-1)+n+1+(j-1)*t_x:n+i*(2*n-1)+(j-1)*t_x)=ones(1,2*n-1); 
   A_eq3(i+(j-1)*n,(2*n-1)*(i-1)+n+1+(j-1)*t_x+(2*n-1)*(n+1-i)+(i-1)*(2*n-2):n+i*(2*n-1)+(j-1)*t_x+(2*n-1)*(n+1-i)+(i-1)*(2*n-2)-1)=-ones(1,2*n-2);
   A_eq3(i+(j-1)*n,j*t_x-n+i)=-1; 
    end
end

A_eq4=sparse(2*n*m,Num_var);
b_eq4=sparse(2*n*m,1);

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


A_ineq1=sparse((4*n*n-n)*m, Num_var);
b_ineq1=sparse((4*n*n-n)*m,1);
A_ineq2=sparse((4*n*n-n)*m, Num_var);
b_ineq2=sparse((4*n*n-n)*m,1);
%%%%Adding another inequality to solve the load matrix inconsistency%%
A_ineqlast=sparse((4*n*n-n)*m, Num_var);
b_ineqlast=sparse((4*n*n-n)*m,1);
%%%%%Loads at depots are zero%%%%
A_eqlast=sparse(2*m,Num_var);
b_eqlast=sparse(2*m,1);
for k=1:m
A_eqlast(k,m*t_x+m*(2*n+2)+(k-1)*(2*n+2)+1)=1;
A_eqlast(k+m,m*t_x+m*(2*n+2)+k*(2*n+2))=1;
end
%%%%%%%%%%


%%%%%%%%


%%%%%%% Setting up big numbers M and W

for k=1:m
           
        for s=1:t_x
          im=arc_set(s,1);
             jm=arc_set(s,2);
                 d_i=d_vec(im+1,1);
                 l_i=l_vec(im+1,1);
                 t_ij=t_vec(s,1);
                 e_j=e_vec(jm+1,1);
                
           M_vec(s+(k-1)*t_x,1)=max(0,l_i+d_i+t_ij-e_j);
       Q_k=Cap_vec(k,1);
       q_i=q_vec(im+1,1);
            W_vec(s+(k-1)*t_x,1)=min(Q_k,Q_k+q_i);    
       
        end
  
end

% M_vec=12000*ones((4*n*n-n)*m,1);
 %W_vec=200*ones((4*n*n-n)*m,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Another approache to calculate A_ineq1 and Aineq2
for cou=1:t_x*m
    
    kp= fix(cou/t_x);
    kpt=rem(cou,t_x);
    if kpt==0
       kpt=t_x; 
       kp=kp-1;
    end
   
    A_ineq1(cou,cou)=M_vec(cou,1);
    
    A_ineq2(cou,cou)=W_vec(cou,1);
    A_ineqlast(cou,cou)=W_vec(cou,1);
    
        ip=arc_set(kpt,1);
    jp=arc_set(kpt,2);
     A_ineq1(cou,t_x*m+ip+1+kp*(2*n+2))=1;
     A_ineq1(cou,t_x*m+jp+1+kp*(2*n+2))=-1;
     
      
      
     A_ineq2(cou,t_x*m+ip+1+kp*(2*n+2)+m*(2*n+2))=1;
     A_ineq2(cou,t_x*m+jp+1+kp*(2*n+2)+m*(2*n+2))=-1;
     
     A_ineqlast(cou,t_x*m+ip+1+kp*(2*n+2)+m*(2*n+2))=-1;
     A_ineqlast(cou,t_x*m+jp+1+kp*(2*n+2)+m*(2*n+2))=1;
     
    
     b_ineq1(cou,1)=M_vec(cou,1)-d_vec(ip+1,1)-t_vec(kpt,1);
     b_ineq2(cou,1)=W_vec(cou,1)-q_vec(jp+1,1);
     
     b_ineqlast(cou,1)=W_vec(cou,1)+q_vec(jp+1,1);
end



     
A_ineq3=sparse(m*n,Num_var);
b_ineq3=sparse(n*m,1);
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



A_ineq4=sparse(m,Num_var);
b_ineq4=T_vec;

for k=1:m
    
 A_ineq4(k,m*t_x+(k-1)*(2*n+2)+1)=-1;   
 A_ineq4(k,m*t_x+k*(2*n+2))=1;   
  
end


%%%%%% Bounds on Variables%%%
lb=sparse(Num_var,1);
ub=inf*ones(Num_var,1);
ub(1:t_x*m,1)=ones(t_x*m,1);




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

A_ineq=sparse([A_ineq1;A_ineq2;A_ineq4]);
b_ineq=[b_ineq1;b_ineq2;b_ineq4];

A_eq=sparse([A_eq1;A_eq2;A_eq3;A_eq4;A_ineq3;A_eqlast]);
b_eq=[b_eq1;b_eq2;b_eq3;b_eq4;b_ineq3;b_eqlast];

% int_var_1=1:m*t_x;
% int_var_2=m*t_x+m*(2*n+2)+1:m*t_x+2*m*(2*n+2);
% int_var=[int_var_1,int_var_2];
%%%%%Using spars structure%%%
% A_eq=sparse(A_eq);
% A_ineq=sparse(A_ineq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%options = optimoptions('intlinprog','Display','iter','OutputFcn',@savemilpsolutions,'PlotFcn',@optimplotmilp);
% options = optimoptions('intlinprog','Display','off','CutMaxIterations',25,'CutGeneration','advanced','IntegerPreprocess','advanced','HeuristicsMaxNodes',100);
% [x,fval,exitflag,output]  = intlinprog(f_vec,int_var,A_ineq,b_ineq,A_eq,b_eq,lb,ub,options)

%%%%Solving with Cplex%%%%
options = cplexoptimset;
   options.Display = 'off';
   options.TolXInteger=1e-2;
   options.mip.tolerances.mipgap=.3;
   options.mip.tolerances.absmipgap=1e-2;
    options.mip.TolFun=.1;
%    options.mip.limits.treememory=15000;
   options.output.intsolfileprefix='on';
   %options.mip.strategy.search=3;
    ctype_vec=zeros(Num_var,1);  
    
   for ictype=1:Num_var
       if (ictype<=m*t_x)
           ctype_vec(ictype,1) = 66;
       elseif ((ictype<=m*t_x+2*m*(2*n+2))&&(ictype>=m*t_x+m*(2*n+2)+1))
           ctype_vec(ictype,1) = 73;
       else
           ctype_vec(ictype,1) = 67;
       end
   end
    
   ctype=char(ctype_vec);

   ctype=ctype';
   %ctypep=ctype(1,1:Num_var);

   %%%%%Freeing some memory%%
   %%clear A_ineq1 A_ineq2 A_ineq3 A_ineq4 A_eq1 A_eq2 A_eq3 A_eq4  b_ineq1 b_ineq2 b_ineq3 b_ineq4 b_eq1 b_eq2 b_eq3 b_eq4 ctype_vec
   
  tic;
   t_cpu = cputime;  
   
 [x,fval,exitflag,output] = cplexmilp(f_vec,A_ineq,b_ineq,A_eq,b_eq,[],[],[],lb,ub,ctype,[ ], options)
  total_time=toc
  e_cpu = cputime-t_cpu
  
%  fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
%    fprintf ('Solution value = %f \n', fval);
%    pause
%%%% using another solver%%%%
% yidx=true(m*t_x,1);
% yidx(m*t_x+1:Num_var,1)=sparse(Num_var-m*t_x,1);
% x=miprog(f_vec,A_ineq,b_ineq,A_eq,b_eq,lb,ub,yidx)

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
        arcs_mat(1:t_x,k)=Passed_arcs((k-1)*t_x+1:k*t_x,1);
        service_start_mat(1:2*n+2,k)=Service_start((k-1)*(2*n+2)+1:k*(2*n+2),1);
        load_mat(1:2*n+2,k)=Loads((k-1)*(2*n+2)+1:k*(2*n+2),1);
        ride_time_mat(1:n,k)=Ride_time(1+(k-1)*n:k*n,1);
    end
    
    %%%%Picking up the relevent results%%%%%
     rel_service_startmat= -ones(2*n+2,m);

    
    for k=1:m
       for i=1:t_x
            if (arcs_mat(i,k)<=1.01 && arcs_mat(i,k)>=.99)
                
                ip=arc_set(i,1)+1;
                jp= arc_set(i,2)+1;

                 rel_service_startmat(ip,k)=service_start_mat(ip,k);
                 rel_service_startmat(jp,k)=service_start_mat(jp,k);

            end
        end
        
    end
   
    %%%%Finding sequence%%%
 
   for k=1:m
       
       [lp,op]=sort(rel_service_startmat(:,k));
       counter=1;
    
       for i=1:2*n+2
           if lp(i,1)>=0
               seq_served_nodes{counter,k}=op(i,1)-1;
               rel_load_mat{counter,k}=load_mat(op(i,1),k);
               rel_service_start_mat{counter,k}=rel_service_startmat(op(i,1),k);
               
               if (op(i,1)<=n+1&&op(i,1)>=2)
                   rel_ride_time_mat{op(i,1)-1,k}=ride_time_mat(op(i,1)-1,k);
                   
                 
               end
               
               counter=counter+1;
           end
       end
   end
   
    seq_served_nodes
    rel_service_start_mat
    rel_load_mat 
    rel_ride_time_mat
    
    %%%%Writing results in Excel%%

%cd('D:/VRP/DARP');
%%%%for linux
%cd('/media/kaveh/DATA/VRP/DARP');
%%%%%%%%%%%%%%%%%%


delete('Results_from_exact.xls') 
header={'Sequence of served nodes','Loads on Vehicles','Starting time','Ride time','Total cost'};
xlswrite('Results_from_exact',header(1),1,'A1');
range_lp = sprintf('A%i', m+2);
    xlswrite('Results_from_exact',header(2),1,range_lp);
    range_lp = sprintf('A%i', 2*m+3);
    xlswrite('Results_from_exact',header(3),1,range_lp);
    range_lp = sprintf('A%i', 3*m+4);
    xlswrite('Results_from_exact',header(4),1,range_lp);
    range_lp = sprintf('A%i', 4*m+5);
    xlswrite('Results_from_exact',header(5),1,range_lp);
    range_lp = sprintf('A%i', 4*m+6);
    xlswrite('Results_from_exact',fval,1,range_lp);

for i=1:m
range_s = sprintf('B%i', i+1);
range_l = sprintf('B%i', i+m+2);
range_b = sprintf('B%i', i+2*m+3);
range_r = sprintf('B%i', i+3*m+4);
    xlswrite('Results_from_exact',seq_served_nodes(:,i)',1,range_s);
    xlswrite('Results_from_exact',rel_load_mat(:,i)',1,range_l);
    xlswrite('Results_from_exact',rel_service_start_mat(:,i)',1,range_b);
    xlswrite('Results_from_exact',rel_ride_time_mat(:,i)',1,range_r);
    
end

        %%%%%%%Visualization
        
        %point_coordinate=10*rand(2,2*n+2);
    point_coordinate=[x_coord y_coord]';
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
                             if (.95<=arcs_mat(i,k)&&arcs_mat(i,k)<=1.111111)
                                plot(xval_mat{i},y_mat{i})
        
       
                               % hold on
                             end
        
        end  
    end
    

end