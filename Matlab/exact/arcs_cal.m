function arcs_set=arcs_cal(n)

arcs_set=zeros(4*n*n-n,2);
counter_arc=0;
for ip=1:(2*n+2)
    
    
    for jp=1:(2*n+2)
        
        if (ip==jp)||(ip==jp+n)||(ip==1&&jp>n+1)|| (jp==2*n+2&&ip<=n+1)||(jp==1&&ip>1)||(ip<n+1&&jp==2*n+2)||(ip==2*n+2)
            
%         elseif jp==2*n+2
%             countp=counter_arc+(n-(ip-n-1))*(2*n-2)+(ip-n-1-1);
%                 arcs_set(countp,1:2)=[ip-1 jp-1];
             
        else
            counter_arc=counter_arc+1;
            arcs_set(counter_arc,1:2)=[ip-1 jp-1];
        end
    end
end

counterp=0;
arcp=zeros(n,2);

for i=1:counter_arc
    
   if arcs_set(i,1)>=n+1 && arcs_set(i,2)==2*n+1 
    counterp=counterp+1;
  
        
     arcs_set(counter_arc+counterp,1:2)=arcs_set(i,1:2);
     arcs_set(i,1:2)=[0 0];
   end
    
end
for i=1:counter_arc
    
   if arcs_set(i,1)==0 && arcs_set(i,2)==0 
        
    arcs_set(i,:)=[];
   
   end
    
end





