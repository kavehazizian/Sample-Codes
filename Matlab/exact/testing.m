function L=LP_TSP(X,Y,Depot)

%[X,Y] are coodinates of nodes distributed in a plane. Depot is the set of
%coordinates of vehicles.
% This program generates feasible routes for vehicles to cover all nodes in
% the plane. The goal is to minimize the max length of a single tour.

[N,N1]=size(X);
[M,M2]=size(Depot);
rr=2*pi*rand(M,1);
City=[X,Y];
drift=[sin(rr) cos(rr)];
Depot=Depot+Depot.*drift/50;

%----------solving LP for clustering -----------
Ca=ceil(N/M)*ones(M-1,1);
Ca(M)=N-sum(Ca);
Z=copl1main1(X,Y,Depot,Ca);

%--------clustering -------
cluster=zeros(ceil(N/M),M);
count=zeros(M,1);
for index=1:N
   i=floor(Z(index,1)/N)+1;
   if (Z(index,1)==(i-1)*N)
       i = i-1;
   end;
   count(i)=count(i)+1;
   cluster(count(i),i)=Z(index,1)-(i-1)*N;
end;

Depot=Depot./(ones(size(Depot))+drift/50);

figure
for i=1:M
   plot(Depot(i,1), Depot(i,2), '*','color','r');
   hold on;
end

for i=1:N
   plot(City(i,1), City(i,2), 'o');
   hold on;
end
%------
%---------------TSP------
L=zeros(1,M);
for i=1:M
    XY=City(cluster(1:count(i),i),:);
    [L(1,i),list]=tsp_con(XY(:,1),XY(:,2), Depot(i,:));
    plotcolor = rand(1,3);
    tempx=[Depot(i,1);XY(:,1)];
    tempy=[Depot(i,2);XY(:,2)];
    plot([tempx(list);Depot(i,1)],[tempy(list);Depot(i,2)],'color',plotcolor);
    hold on;
end




function [tsp_length,list] = tsp_con(X,Y,depot)
        Xtrans=X-depot(1,1);
        Ytrans=Y-depot(1,2);
        [tsp_length,list]= tsp_solve([0;Xtrans],[0;Ytrans],2);
end

%%%%%%%%%%%%BEGIN TSP-METHOD SPECIFIC FUNCTIONS


function [total_path_length,num_cars] = vrp_tsp(Xtrans,Ytrans,k,normchoice)
%solves the simple TSP-based vrp problem, with points located at Xtrans and
%Ytrans and the depot located at the origin, with capacity k and normchoice
%1,2,3 as specified above

%first, get a TSP solution
[path_length,list] = tsp_solve([0;Xtrans],[0;Ytrans],normchoice);
[path_length,list] = tsp_solve(Xtrans,Ytrans,normchoice);
%next, break into component paths
cpt_list = break_into_paths(Xtrans,Ytrans,list,k,normchoice);
%have a list of the startpoints of the path, so plot everything
hold on;
total_path_length = 0;
for i = 1:length(cpt_list)-1
    Xprime = [0;Xtrans(list(cpt_list(i):(cpt_list(i+1)-1)));0];
    Yprime = [0;Ytrans(list(cpt_list(i):(cpt_list(i+1)-1)));0];
    if rem(i,5) == 1
        plot(Xprime,Yprime);
    elseif rem(i,5) == 2
        plot(Xprime,Yprime,'r');
    elseif rem(i,5) == 3
        plot(Xprime,Yprime,'y');
    elseif rem(i,5) == 4
        plot(Xprime,Yprime,'k');
    else
        plot(Xprime,Yprime,'g');
    end
    total_path_length = total_path_length + compute_path_length(Xprime,Yprime);
end

Xprime = [0;Xtrans(list(cpt_list(end):end));0];
Yprime = [0;Ytrans(list(cpt_list(end):end));0];
total_path_length = total_path_length + compute_path_length(Xprime,Yprime);
i = i+1;
    if rem(i,5) == 1
        plot(Xprime,Yprime);
    elseif rem(i,5) == 2
        plot(Xprime,Yprime,'r');
    elseif rem(i,5) == 3
        plot(Xprime,Yprime,'y');
    elseif rem(i,5) == 4
        plot(Xprime,Yprime,'k');
    else
        plot(Xprime,Yprime,'g');
    end
num_cars = i;
title('TSP-based solution');
end





function cpt_list = break_into_paths(Xtrans,Ytrans,list,k,normchoice)
%cpt_list is the starting points of each subpath in 'list' that obeys the
%constraint k; in other words, we have a TSP tour 'list', and a capacity k,
%so this outputs a sequence of paths with length not exceeding k
cpt_list = [];
N = length(Xtrans);
A = distances([[Xtrans;0],[Ytrans;0]]);
endpoint = 0;
IS_FINISHED = false;
while endpoint < length(Xtrans)
    startpoint = endpoint + 1;
    endpoint = startpoint;
    path_dist = A(list(startpoint),N+1);
    while path_dist + A(list(endpoint),N+1) < k & IS_FINISHED == false
        if endpoint == length(Xtrans)
            IS_FINISHED = true;
        else
            endpoint = endpoint + 1;
            path_dist = path_dist + A(list(endpoint-1),list(endpoint));
        end
    end
    if IS_FINISHED == false
    endpoint = endpoint - 1;
    end
    cpt_list = [cpt_list,startpoint];
end

end



%%%%%%%%%%%%%END TSP-METHOD SPECIFIC FUNCTIONS


%%%%%%%%%%%%%BEGIN GENERIC HELPER FUNCTIONS


function path_length = compute_path_length(Xtrans,Ytrans)
%given a set of vectors, computes the path length
A = distances([Xtrans,Ytrans]);
path_length = 0;
for i = 1:length(Xtrans)-1
    path_length = path_length + A(i,i+1);
end
end



function [p,path_length] = lin_kernighan(Xtrans,Ytrans,normchoice)
%computes the lin-kernighan heuristic path length of a TSP tour of Xtrans
%and Ytrans
A = distances([Xtrans,Ytrans]);
matrix_to_tsp(Xtrans,Ytrans,normchoice,'lk.tsp');
!linkern -o output lk.tsp
fid = fopen('output');
C = textscan(fid,'%s');
C = char(C{1});
C = C(:,:);
for i = 1:size(C,1)
    p(i) = str2double(C(i,:));
end
p = p(3:3:end-2)+1;
path_length = 0;
for i = 2:length(p)
    path_length = path_length+A(p(i),p(i-1));
end
path_length = path_length + A(p(1),p(end));
fclose(fid);
!del lk.tsp output
end

function matrix_to_tsp(Xtrans,Ytrans,normchoice,name)
%takes a list of points Xtrans,Ytrans and a choice of norm, and a string
%'name', and makes the appropriate tsp file

n = length(Xtrans);
fid = fopen(name,'w');
fprintf(fid, 'NAME: concorde100\n');
fprintf(fid, 'TYPE: TSP\n');
fprintf(fid, 'DIMENSION: %d\n',n);

if normchoice == 1
    fprintf(fid, 'EDGE_WEIGHT_TYPE: MAN_2D\n');
end

if normchoice == 2
    fprintf(fid, 'EDGE_WEIGHT_TYPE: EUC_2D\n');
end

if normchoice == 3
    fprintf(fid, 'EDGE_WEIGHT_TYPE: MAX_2D\n');
end

fprintf(fid, 'NODE_COORD_SECTION\n');
for i = 1:n
    fprintf(fid, '%d %f %f\n',i,Xtrans(i),Ytrans(i));
end
fclose(fid);
end




function [path_length,list] = tsp_solve(Xtrans,Ytrans,normchoice)
%Takes some vectors Xtrans, Ytrans, and a norm choice, and solves the TSP
%problem with them
A = distances([Xtrans,Ytrans]); %eventually, this needs to be updated to reflect different norms
matrix_to_tsp(Xtrans,Ytrans,2,'tosolve.tsp');
!concorde tosolve.tsp
list=sol_to_matrix('tosolve.sol');

%need the length of the path
path_length = 0;
for i = 2:length(list)
    path_length = path_length + A(list(i),list(i-1));
end
path_length = path_length + A(list(1),list(end));
!del tosolve.tsp Otosolve.pul Otosolve.sav Otosolve.mas tosolve.mas tosolve.pul tosolve.sav tosolve.sol
end




function list = sol_to_matrix(name)
%makes a list of numbers for a .sol file, to be used by concorde
fid = fopen(name);
C = textscan(fid,'%s');
C = char(C{1});
C = C(2:end,:);
for i = 1:size(C,1)
    list(i) = str2double(C(i,:));
end
fclose(fid);
list=list+1;
end
end