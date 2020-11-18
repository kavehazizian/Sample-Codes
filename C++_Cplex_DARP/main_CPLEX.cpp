#include <iostream>
#include <sstream>
#include <iomanip>
#include<fstream>
#include<string>
#include <vector>
#include <armadillo>
#include <algorithm>
#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <limits>
#include <time.h>
#include </media/kaveh/DATA/Installed_cplex_linux/cplex/include/ilcplex/ilocplex.h>
//#include </media/kaveh/DATA/Installed_cplex_linux/concert/ilconcert/ilomodel.h>
//#include </media/kaveh/DATA/Installed_cplex_linux/concert/ilconcert/iloalg.h>
//#include "/media/kaveh/DATA/C++_Cplex_DARP/Reading_inputs/StanfordCPPLib/console.h"

using namespace std;




arma::umat arcs_cal(int n) {
    arma::umat arcs_set;
    arcs_set.zeros(4*n*n,2);
    int counter_arc=0;


    for (int ip = 1; ip <= 2*n+2; ip++) {

    for (int jp = 1; jp <= 2*n+2; jp++) {
    if ( (ip==jp)||(ip==jp+n)||(ip==1&&jp>n+1)|| (jp==2*n+2&&ip<=n+1)||(jp==1&&ip>1)||(ip<n+1&&jp==2*n+2)||(ip==2*n+2) )

{}

 else
            {
            arcs_set(counter_arc,0)=ip-1;
             arcs_set(counter_arc,1)=jp-1;

counter_arc=counter_arc+1;
             }

    }
    }

//cout<<counter_arc<<endl;
 int counterp=0;

for (int i=0;i<counter_arc;i++)
{
   if (arcs_set(i,0)>=n+1 && arcs_set(i,1)==2*n+1)
    {


     arcs_set(counter_arc+counterp,arma::span(0,1))=arcs_set(i,arma::span(0,1));
     arcs_set(i,0)= 0;
     arcs_set(i,1)= 0;
     //cout<<i<<endl;
     counterp=counterp+1;

   }

}

for (int i=0;i<counter_arc;i++){

   if (arcs_set(i,0)==0 && arcs_set(i,1)==0)

    {arcs_set.shed_row(i);

   }

}


//arma::mat arcp.zeros(n,2);




    return arcs_set;

}









int main()
{

struct winsize w;

std::cout.precision(2);



ifstream myfile("darp_a_2_20.txt");
int m=2,n=20;
double Max_ride=30;
int i=0;
int j,t_x;
int no;
double x_c,y_c,d_v,q_v,e_v,l_v;

t_x=4*n*n-n;
arma::vec nodes(2*n+2),x_coord(2*n+2),y_coord(2*n+2),d_vec(2*n+2), q_vec(2*n+2),e_vec(2*n+2), l_vec(2*n+2),Cap_vec(m),T_vec(m),t_vec(t_x),cost_matrix_vec;
t_vec.zeros(t_x);
cost_matrix_vec.zeros(t_x*m);


Cap_vec=3*Cap_vec.ones(m);
T_vec=600*T_vec.ones(m);

while (myfile>>no>>x_c>>y_c>>d_v>>q_v>>e_v>>l_v){


x_coord(no)=x_c;
y_coord(no)=y_c;
d_vec(no)=d_v;
q_vec(no)=q_v;
e_vec(no)=e_v;
l_vec(no)=l_v;


}


arma::umat arc_set(4*n*n-n,2);

arc_set=arcs_cal(n);

int ip,jp;
for (int i=0;i<t_x;i++){
     ip=arc_set(i,0);
     jp=arc_set(i,1);
    t_vec(i)= sqrt(pow(x_coord(ip)-x_coord(jp),2)+pow(y_coord(ip)-y_coord(jp),2));
}
for (int k=1;k<=m;k++){
    cost_matrix_vec.rows((k-1)*t_x,k*t_x-1)=t_vec;

  }


 /////Matrices and Vectors constructions start here///
int Num_var=4*m*(n*n+n+1);







/// objective fun

arma:: vec f_vec;
f_vec.zeros(Num_var);
f_vec.rows(0,(4*n*n-n)*m-1)=cost_matrix_vec;

///Equality constraints
arma::sp_mat A_eq1(n,Num_var);

arma::vec b_eq1(n);
b_eq1.ones();



for (int i=1;i<=n;i++){
    for (int j=1;j<=m;j++){
    A_eq1(i-1,arma::span((2*n-1)*(i-1)+n+1+(j-1)*t_x-1,n+i*(2*n-1)+(j-1)*t_x-1)).ones();}
    }

arma::sp_mat A_eq2(2*m,Num_var);
arma::vec b_eq2(2*m);
b_eq2.ones();
////starting elements for delivery arcs
int t_n=4*n*n-2*n;
for (int i=1;i<=m;i++){
A_eq2(i-1,arma::span(t_x*(i-1),t_x*(i-1)+n-1)).ones();
A_eq2(m+i-1,arma::span(t_x*(i-1)+t_n,t_x*(i-1)+t_n+n-1)).ones();
}

arma::sp_mat A_eq3(n*m,Num_var);
arma::vec b_eq3(n*m),A_eq3_vec(2*n-2);
b_eq3.zeros();
A_eq3_vec.ones();
for (int j=1;j<=m;j++){
    for (int i=1;i<=n;i++){

   A_eq3(i+(j-1)*n-1,arma::span( (2*n-1)*(i-1)+n+1+(j-1)*t_x-1,n+i*(2*n-1)+(j-1)*t_x-1)).ones();
   A_eq3(i+(j-1)*n-1,arma::span( (2*n-1)*(i-1)+n+1+(j-1)*t_x+(2*n-1)*(n+1-i)+(i-1)*(2*n-2)-1,n+i*(2*n-1)+(j-1)*t_x+(2*n-1)*(n+1-i)+(i-1)*(2*n-2)-2))=-trans(A_eq3_vec);
   A_eq3(i+(j-1)*n-1,j*t_x-n+i-1)=-1;
    }
}


arma::sp_mat A_eq4(2*n*m,Num_var);
arma:: vec b_eq4(2*n*m);
b_eq4.zeros();


for (int k=1;k<=m;k++){
    for (int i=1;i<=2*n;i++){

        for (int s=1;s<=t_x;s++){
        if   (arc_set(s-1,0)==i){

            A_eq4(i+(k-1)*2*n-1,s+(k-1)*t_x-1)=1;}
        else if (arc_set(s-1,1)==i){
            A_eq4(i+(k-1)*2*n-1,s+(k-1)*t_x-1)=-1;
        }
        }
    }
}

arma::sp_mat A_ineq1((4*n*n-n)*m, Num_var),A_ineq2((4*n*n-n)*m, Num_var),A_ineqlast((4*n*n-n)*m, Num_var),A_eqlast(2*m,Num_var);
arma::vec b_ineq1((4*n*n-n)*m),b_ineq2((4*n*n-n)*m),b_ineqlast((4*n*n-n)*m),b_eqlast(2*m),M_vec((4*n*n-n)*m),W_vec((4*n*n-n)*m);

//Adding another inequality to solve the load matrix inconsistency%%

//Loads at depots are zero%%%%
b_ineq1.zeros();
b_ineq2.zeros();
b_ineqlast.zeros();
b_eqlast.zeros();
for (int k=1;k<=m;k++){
A_eqlast(k-1,m*t_x+m*(2*n+2)+(k-1)*(2*n+2))=1;
A_eqlast(k+m-1,m*t_x+m*(2*n+2)+k*(2*n+2)-1)=1;
}

// Setting up big numbers M and W
M_vec.ones();
W_vec.ones();
for (int k=1;k<=m;k++){

        for (int s=1;s<=t_x;s++){
         int            im=arc_set(s-1,0);
         int            jm=arc_set(s-1,1);
         double        d_i=d_vec(im);
         double        l_i=l_vec(im);
         double       t_ij=t_vec(s-1);
         double        e_j=e_vec(jm);

           M_vec(s+(k-1)*t_x-1)=max<double>(0,l_i+d_i+t_ij-e_j);
       double Q_k=Cap_vec(k-1);
       double q_i=q_vec(im);
            W_vec(s+(k-1)*t_x-1)=min<double>(Q_k,Q_k+q_i);

        }

}


//// Another approache to calculate A_ineq1 and Aineq2
for (int cou=1;cou<=t_x*m;cou++){

    int kp= (int)cou/t_x;
    int kpt=cou%t_x;
    if (kpt==0){
    kpt=t_x;
    }


    A_ineq1(cou-1,cou-1)=M_vec(cou-1);

    A_ineq2(cou-1,cou-1)=W_vec(cou-1);
    A_ineqlast(cou-1,cou-1)=W_vec(cou-1);
    if (kp==m){
        kp=m-1;
    }
        ip=arc_set(kpt-1,0);
    jp=arc_set(kpt-1,1);
     A_ineq1(cou-1,t_x*m+ip+kp*(2*n+2))=1;
     A_ineq1(cou-1,t_x*m+jp+kp*(2*n+2))=-1;



     A_ineq2(cou-1,t_x*m+ip+kp*(2*n+2)+m*(2*n+2))=1;
     A_ineq2(cou-1,t_x*m+jp+kp*(2*n+2)+m*(2*n+2))=-1;

     A_ineqlast(cou-1,t_x*m+ip+kp*(2*n+2)+m*(2*n+2))=-1;
     A_ineqlast(cou-1,t_x*m+jp+kp*(2*n+2)+m*(2*n+2))=1;


     b_ineq1(cou-1)=M_vec(cou-1)-d_vec(ip)-t_vec(kpt-1);
     b_ineq2(cou-1)=W_vec(cou-1)-q_vec(jp);

     b_ineqlast(cou-1)=W_vec(cou-1)+q_vec(jp);
}

arma::sp_mat A_ineq3(m*n,Num_var);
arma:: vec b_ineq3(n*m,1);
b_ineq3.zeros();
int Num_ineq3=m*(4*n*n+3*n+4);

int cout_3=0;
for (int k=1;k<=m;k++){

    for (int i=1;i<=n;i++){

 A_ineq3(cout_3,Num_ineq3+i+(k-1)*n-1)=-1;
  A_ineq3(cout_3,m*t_x+i+n+1+(k-1)*(2*n+2)-1)=1;
  A_ineq3(cout_3,m*t_x+i+1+(k-1)*(2*n+2)-1)=-1;
  b_ineq3(cout_3)=d_vec(i);
 cout_3=cout_3+1;
    }
}

arma::sp_mat A_ineq4(m,Num_var);
arma::vec b_ineq4(m);
b_ineq4=T_vec;

for (int k=1;k<=m;k++){

 A_ineq4(k-1,m*t_x+(k-1)*(2*n+2))=-1;
 A_ineq4(k-1,m*t_x+k*(2*n+2)-1)=1;

}

///// Bounds on Variables%%%
arma::vec lb(Num_var),ub(Num_var),one_vec(m*n);
one_vec.ones();
lb.zeros();
double a_inf = std::numeric_limits<double>::infinity();
ub=a_inf*ub.ones();
ub(arma::span(0,t_x*m-1)).ones();


for (int k=1;k<=m;k++){
    lb(arma::span (t_x*m+(k-1)*(2*n+2),t_x*m+k*(2*n+2)-1))=e_vec;
    ub(arma::span (t_x*m+(k-1)*(2*n+2),t_x*m+k*(2*n+2)-1))=l_vec;
    for (int i=1;i<=n;i++){
       lb(Num_var-m*n+i+(k-1)*n-1)=t_vec(n+(i-1)*(2*n-1)+(n-1)+i-1);

    }

    ub(arma::span(Num_var-m*n,Num_var-1)) =Max_ride*one_vec;

    for (int ip=1;ip<=2*n+2;ip++){
        lb(Num_var-m*n-m*(2*n+2)+(k-1)*(2*n+2)+ip-1)=max<double>(0,q_vec(ip-1));
        ub(Num_var-m*n-m*(2*n+2)+(k-1)*(2*n+2)+ip-1)=min<double>(Cap_vec(k-1),Cap_vec(k-1)+q_vec(ip-1));

    }



}
//// Assembling Matrices
arma::sp_mat A_ineq,A_ineqp,A_ineqpp, A_eq,A_eqp,A_eqpp,A_eqppp,A_eqf;
arma::vec b_ineq,b_ineqp,b_ineqpp,b_eq,b_eqp,b_eqpp,b_eqppp,b_eqf;
 	A_ineqp= join_vert(A_ineq1,A_ineq2);
    A_ineqpp=join_vert(A_ineq4,A_ineqlast);
    A_ineq=join_vert(A_ineqp,A_ineqpp);

b_ineqp=join_vert(b_ineq1,b_ineq2);
b_ineqpp=join_vert(b_ineq4,b_ineqlast);
b_ineq=join_vert(b_ineqp,b_ineqpp);

    A_eqp= join_vert(A_eq1,A_eq2);
    A_eqpp=join_vert(A_eq3,A_eq4);
    A_eqppp=join_vert(A_ineq3,A_eqlast);
    A_eqf=join_vert(A_eqp,A_eqpp);
    A_eq=join_vert(A_eqf,A_eqppp);


b_eqp= join_vert(b_eq1,b_eq2);
b_eqpp=join_vert(b_eq3,b_eq4);
b_eqppp=join_vert(b_ineq3,b_eqlast);
b_eqf=join_vert(b_eqp,b_eqpp);
b_eq=join_vert(b_eqf,b_eqppp);
/////Freeing memory////////////////////////
A_ineq1.clear();
A_ineq2.clear();
A_ineq3.clear();
A_ineq4.clear();
A_ineqlast.clear();
A_ineqp.clear();
A_ineqpp.clear();
A_eq1.clear();
A_eq2.clear();
A_eq3.clear();
A_eq4.clear();
A_eqlast.clear();
A_eqf.clear();
A_eqp.clear();
A_eqpp.clear();
A_eqppp.clear();

b_ineq1.clear();
b_ineq2.clear();
b_ineq3.clear();
b_ineq4.clear();
b_ineqlast.clear();
b_ineqp.clear();
b_ineqpp.clear();
b_eq1.clear();
b_eq2.clear();
b_eq3.clear();
b_eq4.clear();
b_eqlast.clear();
b_eqf.clear();
b_eqp.clear();
b_eqpp.clear();
b_eqppp.clear();
//////////////////////////////////

////////CPLEX Starts///////
ILOSTLBEGIN
// creating an environment
	IloEnv env;

	// creating a model
	IloModel mymodel(env);

	// variable declaration
	  //IloNumVar x = cplex.numVar(lb, ub, IloNumVarType.Float);

int c_U,c_W,c_R;
c_U=m*(2*n+2);
c_W=m*(2*n+2);
c_R=n*m;
IloNumVarArray X(env);//;,m*t_x,0,1,ILOINT);
IloNumVarArray U(env);//,c_U,0,IloInfinity,ILOFLOAT);
IloNumVarArray W(env);
IloNumVarArray R(env);

//IloNumVarArray x(env,m*(4*n*n-n));
for (int i=0;i<m*t_x;i++){

X.add(IloNumVar(env,0, 1,ILOINT));

}



for (int i=m*t_x;i<m*t_x+c_U;i++){

U.add(IloNumVar(env,lb(i), ub(i),ILOFLOAT) );

}

for (int i=m*t_x+c_U;i<m*t_x+c_U+c_W;i++){

W.add(IloNumVar(env,lb(i), ub(i),ILOINT)) ;


}

for (int i=m*t_x+c_U+c_W;i<m*t_x+c_U+c_W+c_R;i++){

R.add(IloNumVar(env,lb(i), ub(i),ILOFLOAT));


}
IloExpr expr(env);
for (int i = 0; i <m*t_x; i++){
expr=expr+ X[i]*cost_matrix_vec(i);
//if(i==3){
//cout<<expr<<endl;}

}
mymodel.add(IloMinimize (env,expr));

for(IloInt i=0;i<A_eq.n_rows;i++)
{
IloExpr ctr(env);
for (IloInt j=0;j<A_eq.n_cols;j++)
{ if(j<m*t_x){
	ctr+=A_eq(i,j)*X[j];}
	else if(j<m*t_x+c_U&& j>=m*t_x){ctr+=A_eq(i,j)*U[j-m*t_x];}
	else if(j>=m*t_x+c_U&& j<m*t_x+c_U+c_W){ctr+=A_eq(i,j)*W[j-m*t_x-c_U];}
	else{ctr+=A_eq(i,j)*R[j-m*t_x-c_U-c_W];}
	}
mymodel.add(ctr==b_eq(i));
}

for(IloInt i=0;i<A_ineq.n_rows;i++)
{
IloExpr ctr_ineq(env);
for (IloInt j=0;j<A_ineq.n_cols;j++)
{ if(j<m*t_x){
	ctr_ineq+=A_ineq(i,j)*X[j];}
	else if(j<m*t_x+c_U&& j>=m*t_x){ctr_ineq+=A_ineq(i,j)*U[j-m*t_x];}
	else if(j>=m*t_x+c_U&& j<m*t_x+c_U+c_W){ctr_ineq+=A_ineq(i,j)*W[j-m*t_x-c_U];}
	else{ctr_ineq+=A_ineq(i,j)*R[j-m*t_x-c_U-c_W];}
	}
mymodel.add(ctr_ineq<=b_ineq(i));
}
IloCplex mycplex(mymodel);

// solving the problem
 clock_t t_start, t_end;
 t_start = clock();
 mycplex.setOut(env.getNullStream());
mycplex.solve();
mycplex.out() << "Solution Status is " << mycplex.getStatus() << endl;
t_end = clock() ;
cout<<"Total_time="<<(t_end-t_start)/ CLOCKS_PER_SEC<<"Second\n \n";
cout<< " Total cost is "<<fixed<< mycplex.getObjValue() << "\n\n\n";

mycplex.setWarning(env.getNullStream());
mycplex.setError(env.getNullStream());


mycplex.setParam(IloCplex::MIPDisplay, 0);//MIP node log display information-No display until optimal solution has been found
mycplex.setParam(IloCplex::SimDisplay, 0);//No iteration messages until solution
mycplex.setParam(IloCplex::BarDisplay, 0);//No progress information
mycplex.setParam(IloCplex::NetDisplay, 0);//Network logging display indicator
mycplex.setParam(IloCplex::BrDir, 0);

mycplex.setParam(IloCplex::ParallelMode,0);
if ( mycplex.getStatus() == IloAlgorithm::Infeasible ||
           mycplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded ) {
           env.out() << endl << "*** Model is infeasible ***" << endl;
           ofstream DARP_results;

  DARP_results.open ("DARP_results.txt");
  DARP_results << "Problem is infeasible or Unbounded\n\n";
           }
else{
mycplex.writeSolution("sol.txt");
// the results
IloNumArray val_x(env),val_u(env),val_w(env),val_r(env);
 mycplex.getValues(val_x, X);
 mycplex.getValues(val_u, U);
 mycplex.getValues(val_w, W);
 mycplex.getValues(val_r, R);






//for (int i=0;i<c_R;i++){
//cout<<"x="<<mycplex.getValue(R[i]);

//}


//cout<<mycplex.getValues(X[])<<endl;
/////////////////////////////
//cout<<	Num_var<<endl;

 ////////////////////////////
 // releasing the memory

/////// Visualizing the REsults///
arma::fmat Start_time_mat(2*n+2,m);
Start_time_mat.ones();
Start_time_mat=-Start_time_mat;

int i_passed, j_passed;
for (int k=0;k<m;k++){

for(int i=0;i<t_x;i++){
if(val_x[k*t_x+i]>.99&&val_x[k*t_x+i]<1.01){
i_passed=arc_set(i,0);
j_passed=arc_set(i,1);
Start_time_mat(i_passed,k)=val_u[k*(2*n+2)+i_passed];
Start_time_mat(j_passed,k)=val_u[k*(2*n+2)+j_passed];
//Load_mat(i_passed,k)=val_w[k*(2*n+2)+i_passed];
//Load_mat(j_passed,k)=val_w[k*(2*n+2)+j_passed];
//ridetime_mat(j_passed,k)=val_w[k*n+j_passed]
}

}
}
arma::fvec lp(2*n+2);
arma::uvec op(2*n+2);
arma::Mat<int> seq_served_nodes(m,2*n+2),rel_load_mat(m,2*n+2);
arma::fmat rel_service_start_mat(m,2*n+2),rel_ride_time_mat(m,n);

seq_served_nodes.zeros();
rel_load_mat.zeros();
rel_service_start_mat.zeros();
rel_ride_time_mat.zeros();


//cout<< " Total U "<< val_u<< "\n";
for (int k=0;k<m;k++){

       lp=arma::sort(Start_time_mat.col(k));
       op=sort_index(Start_time_mat.col(k));
       int counter=0;

for (int i=0;i<2*n+2;i++){
           if (lp(i)>=0){
               seq_served_nodes(k,counter)=op(i);

               rel_load_mat(k,counter)=val_w[k*(2*n+2)+op(i)];
               rel_service_start_mat(k,counter)=val_u[k*(2*n+2)+op(i)];


               if (op(i)<n+1&&op(i)>0){

                  rel_ride_time_mat(k,op(i)-1)=val_r[k*n+op(i)-1];
//

               }

               counter=counter+1;
           }
       }
   }
cout.setf(ios::fixed);
//cout << setprecision(0)<< rel_service_start_mat<< "\n";
//cout<< " Total cost is "<< fixed<<rel_service_start_mat<< "\n";





//
ofstream DARP_results;
int iw;
  DARP_results.open ("DARP_results.txt");
  DARP_results  << fixed << noshowpoint;
DARP_results  << setprecision(2);
DARP_results << " Total cost is "<< mycplex.getObjValue() << "\n\n\n";

  for (int k=0;k<m;k++){
  iw=1;

  while (seq_served_nodes(k,iw)!=0 && iw<2*n+2){
iw=iw+1;

}


DARP_results << "Vehicle "<<k+1<<"\n\n";
DARP_results  << "Served nodes:\n"<<seq_served_nodes( arma::span(k,k),arma::span(0,iw-1) )<<"\n";
DARP_results << "Service Starting time at nodes:\n"<<rel_service_start_mat( arma::span(k,k),arma::span(0,iw-1) )<<"\n";
DARP_results << "Loads of Vehicle after leaving nodes:\n"<<rel_load_mat( arma::span(k,k),arma::span(0,iw-1) )<<"\n";
DARP_results <<  "Ride time of the served nodes\n"<<rel_ride_time_mat.row(k)<<"\n \n \n \n";


  }


  DARP_results.close();

//cout<<rel_ride_time_mat<<endl;
}
  env .end ();

    return 0;
}
