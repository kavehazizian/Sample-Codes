function vec=vec_change(mymatrix)
[m,n]=size(mymatrix);
for i=1:m
    vec(n*(i-1)+1:n*i)=mymatrix(i,1:n);
end
vec=vec';