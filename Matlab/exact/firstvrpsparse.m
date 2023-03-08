clc;
clear all;
close all;
tic;
cost_matrix=zeros(6,6);
cost_matrix(1,2:6)=[28,39,50,57,83];
cost_matrix(2,3:6)=[52,49,74,72];
cost_matrix(3,4:6)=[48,48,111];
cost_matrix(4,5:6)=[35,85];
cost_matrix(5,6)=105;
Num_vhiecles=10;
cap=15000;
% creating symmetric cost matrix
[mc,nc]=size(cost_matrix);
%%%%%%Number of variables 
Num_var=2*mc*nc*Num_vhiecles+mc*Num_vhiecles;
f_vec=sparse(Num_var,1);
for i=1:mc
    for j=1:nc
        cost_matrix(j,i)=cost_matrix(i,j);
    end
end
vec_cost=vec_change(cost_matrix);
for i_f=1:Num_vhiecles
    f_vec(mc*nc*(i_f-1)+1:mc*nc*i_f,1)=vec_cost;
end
%%%% Demands
Dem=0.1*[18000,26000,11000,30000,21000]';

%%%%%Inequality constraints%%%%%%%%%%
A_ineq=sparse(nc-1+mc*nc*Num_vhiecles,Num_var);
%%%%Assignment
for kyk=1:Num_vhiecles
for iyk=1:nc-1
    A_ineq(iyk,nc*nc*Num_vhiecles+1+(kyk-1)*nc+iyk)=-1;
end
end
%%%% Vhiecle Capacity
A_ineq(nc-1+1:nc-1+nc*nc*Num_vhiecles,1:2*nc*nc*Num_vhiecles+nc*Num_vhiecles)=[-cap*eye(nc*nc*Num_vhiecles),zeros(nc*nc*Num_vhiecles,nc*Num_vhiecles),eye(nc*nc*Num_vhiecles)];

b_ineq=sparse(nc-1+mc*nc*Num_vhiecles,1);
b_ineq(1:nc-1,1)=-ones(nc-1,1);

%%%%% Equality Constraints%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_eq=sparse(2*mc*Num_vhiecles+(mc-1)*Num_vhiecles+Num_vhiecles,Num_var);
b_eq=sparse(2*mc*Num_vhiecles+(mc-1)*Num_vhiecles+Num_vhiecles,1);
for kieq=1:Num_vhiecles
for ipaeq=1:mc
   %%%%% Leave node
   A_eq(ipaeq+(kieq-1)*mc,(ipaeq-1)*nc+1+(kieq-1)*nc*mc:ipaeq*nc+(kieq-1)*nc*mc)=ones(1,nc); 
   A_eq(ipaeq+(kieq-1)*mc,mc*nc*Num_vhiecles+(kieq-1)*nc+ipaeq)=-1;
   %%%%%Enter node
   A_eq(ipaeq+(kieq-1)*mc+mc*Num_vhiecles,mc*nc*Num_vhiecles+(kieq-1)*nc+ipaeq)=-1;
   A_eq((kieq-1)*mc+mc*Num_vhiecles+1:(kieq-1)*mc+mc*Num_vhiecles+mc,(ipaeq-1)*nc+(kieq-1)*mc+1:ipaeq*nc+(kieq-1)*mc)=eye(mc,nc);
end

%%%% Depot
A_eq(2*mc*Num_vhiecles+(mc-1)*Num_vhiecles+kieq,mc*nc*Num_vhiecles+(kieq-1)*nc+1)=1;


end
%%%% Depot

b_eq(2*mc*Num_vhiecles+(mc-1)*Num_vhiecles+1:2*mc*Num_vhiecles+(mc-1)*Num_vhiecles+Num_vhiecles,1)=ones(Num_vhiecles,1);
%%%% Flow Balance
for kb=1:Num_vhiecles
    beta=mc*nc*Num_vhiecles+mc*Num_vhiecles+(kb-1)*mc*nc;
    betay=mc*nc*Num_vhiecles+(kb-1)*nc;
for ib=2:nc
    alpha=2*mc*Num_vhiecles+(kb-1)*(nc-1)+ib-1;
   
    for jb=1:nc
   A_eq(alpha,beta+(jb-1)*nc+ib)=1;
   if ib==jb
      A_eq(alpha,beta+nc*(ib-1)+1:beta+nc*(ib-1)+nc)=-ones(1,nc);
      A_eq(alpha,beta+nc*(ib-1)+jb)=0;
      A_eq(alpha,betay+ib)=-Dem(ib-1);
   end
   
    
    
    end
end
end

%%%Lower and Upper bonds on variable
lb=zeros(Num_var,1);
ub=Inf*ones(Num_var,1);
%%%%%To ensure that Y and Z are binary variables
ub(1:Num_var-nc*mc*Num_vhiecles,1)=ones(Num_var-nc*mc*Num_vhiecles,1);
%%%%%% To ensure that all variables are binary variables
%%%%% ub=ones(Numvar,1);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Specifying which variables are integers
%intcon=1:Num_var-nc*mc*Num_vhiecles;
intcon=1:Num_var;
%%%%% Solving the MILP%%%%%
[x,fval,exitflag,output] = intlinprog(f_vec,intcon,A_ineq,b_ineq,A_eq,b_eq,lb,ub);
toc