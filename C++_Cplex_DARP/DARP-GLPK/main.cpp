#include <iostream>
#include <sstream>
#include<fstream>
#include <iomanip>
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
#include <lemon/lp.h>
#include <glpk.h>



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



///// g++ -O3 -o main main.cpp -lemon -lglpk


using namespace lemon;




int main()
{

struct winsize w;

std::cout.precision(2);



ifstream myfile("input_darp.txt");
int m=2,n=16;
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
T_vec=480*T_vec.ones(m);
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
    kp=kp-1;
    }


    A_ineq1(cou-1,cou-1)=M_vec(cou-1);

    A_ineq2(cou-1,cou-1)=W_vec(cou-1);
    A_ineqlast(cou-1,cou-1)=W_vec(cou-1);

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
//A_ineq1.clear();
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

//b_ineq1.clear();
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

Mip  mip;
 //arma::vec x(Num_var);

std::vector<Mip::Col> x(Num_var);
for(int i=0;i<Num_var;i++){
  x[i] = mip.addCol();

}

  // Add two columns (variables) to the problem

//mip.addColSet(x);
// Add rows (constraints) to the problem


int row_ineq=A_ineq.n_rows;
// /*
   for(int i=0;i<row_ineq;i++){
Mip::Expr e_ineq;
   for(int j=0;j<Num_var;j++)
 {
 e_ineq += A_ineq(i,j) * x[j];
//mip.addRow(A_ineq(i,j)*x[j] <=b_ineq(i));

  }
  mip.addRow(e_ineq <= b_ineq(i));
  }
// */
int row_eq=A_eq.n_rows;


   for(int i=0;i<row_eq;i++){
    Mip::Expr e_eq;
   for(int j=0;j<Num_var;j++)
 {

   e_eq += A_eq(i,j) * x[j];
  //mip.addRow(A_eq(i,j)*x[j]=b_eq(i));

  }
  mip.addRow(e_eq==b_eq(i));
  }

  // Set lower and upper bounds for the columns (variables)
  for (int i=0;i<Num_var;i++){
  mip.colLowerBound(x[i], lb(i));
  mip.colUpperBound(x[i], ub(i));

  }

  // Set the type of the columns
  for(int i=0;i<m*t_x;i++){
  mip.colType(x[i], Mip::INTEGER);
  }
  for(int i=m*t_x;i<m*(t_x+2*n+2);i++){
  mip.colType(x[i], Mip::REAL);
  }
for(int i=m*(t_x+2*n+2);i<m*(t_x+2*(2*n+2) );i++){
  mip.colType(x[i], Mip::INTEGER);
  }
for(int i=m*(t_x+2*(2*n+2));i<m*(t_x+2*(2*n+2)+n);i++){
  mip.colType(x[i], Mip::REAL);
  }

  // Specify the objective function
 mip.min();
for(int i=0;i<Num_var;i++){
   //mip.obj(mip.obj()+f_vec(i)*x[i]);
  mip.objCoeff(x[i],f_vec(i));
}


// Solve the problem using the underlying MIP solver
  mip.solve();


// Print the results

    cout << "Objective function value: " << mip.solValue() << endl;



switch (mip.type()) {
        case 0: cout << "Feasible solution hasn't been found (but may exist)."<< endl;
         break;
        case 1: cout << "The problem has no feasible solution"<< endl;   //execution starts at this case label
         break;
        case 2: cout << "Feasible solution found."<< endl;
         break;
        case 3: cout <<"Optimal solution exists and found"<< endl;
         break;
        case 4: cout << "The cost function is unbounded. The Mip or at least the relaxed problem is unbounded."<< endl;


    }
//cout<< mip.sol(x[m*t_x-1])<<"\t"<< mip.sol(x[m*t_x])<<"\t"<<mip.sol(x[m*(t_x+2*n+2)-1])<<"\t"<<mip.sol(x[m*(t_x+2*n+2)])<<endl;

  /////////////////////////Printing out results///
  if ( mip.type() == 2 || mip.type() == 3 ) {

           ofstream DARP_results_GLPK;

  DARP_results_GLPK.open ("DARP_results_GLPK.txt");
  DARP_results_GLPK << "Solution found!\n\n";
  DARP_results_GLPK <<  " Total cost is "<< mip.solValue() << "\n\n\n";



arma::fmat Start_time_mat(2*n+2,m);
Start_time_mat.ones();
Start_time_mat=-Start_time_mat;
arma::uvec x_fin(m*t_x), val_w(m*(2*n+2));
arma::fvec val_u(m*(2*n+2)),val_r(m*n);
arma::sp_mat x_test(Num_var,1);
arma::vec x_test_mat(Num_var);




for(int i=0;i<m*t_x;i++){
x_fin[i]=mip.sol(x[i]);

}

for(int i=m*t_x;i<m*(t_x+2*n+2);i++)
{val_u(i-m*t_x)=mip.sol(x[i]);

}
////testing Ineqs
//for(int i=0;i<Num_var;i++){
//x_test(i,0)=mip.sol(x[i]);
//x_test_mat(i)=mip.sol(x[i]);
//}
//int row_ineq1=A_ineq1.n_rows;
//arma::vec test_vec(row_ineq1);
//test_vec=A_ineq1*x_test-b_ineq1;
//for(int i=0;i<row_ineq1;i++){
//if(i==1007){
//cout<<"A_ineq1*x="<<A_ineq1(i,arma::span(0,Num_var-1))*x_test<<endl;
//cout<<"Infaesible Inequlities="<< test_vec(i)<<"\t b_ineq(i)="<<b_ineq(i)<<endl;
//cout<<"Arc="<<arc_set(i,0)<<arc_set(i,1)<<endl;}
////cout<<val_u(arc_set(i,0))<<val_u(arc_set(i,1))<<endl;
//}
//cout<<"valu="<<val_u(2*n+1)<<"\t"<<val_u(2*n)<<endl;
//
//cout<<"tij="<<t_vec(1007)<<endl;
//
/////////////////Transfering data to Matlab///
//ofstream abcdef;
//    abcdef.open("/media/kaveh/DATA/matlab_linux/bin/mat_test.bin",ios::out | ios::trunc | ios::binary);
//    arma::mat A_ineq_test(A_ineq.n_rows,Num_var);
//
//    A_ineq_test=A_ineq;
//
//abcdef<<A_ineq_test<<endl;
////abcdef<<x_test_mat;
//abcdef.close();

/////////////////////////////////
for(int i=m*(t_x+2*n+2);i<m*(t_x+2*(2*n+2) );i++)
{val_w(i-m*(t_x+2*n+2))=mip.sol(x[i]);
}

for(int i=m*(t_x+2*(2*n+2) );i<m*(t_x+2*(2*n+2)+n );i++)
{val_r(i-m*(t_x+2*(2*n+2)) )=mip.sol(x[i]);
}


int i_passed, j_passed;
for (int k=0;k<m;k++){

for(int i=0;i<t_x;i++){
if(x_fin[k*t_x+i]>.99&&x_fin[k*t_x+i]<1.01){
i_passed=arc_set(i,0);
j_passed=arc_set(i,1);
Start_time_mat(i_passed,k)=val_u[k*(2*n+2)+i_passed];
Start_time_mat(j_passed,k)=val_u[k*(2*n+2)+j_passed];

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



   int iw;

  DARP_results_GLPK  << fixed << noshowpoint;
DARP_results_GLPK  << std::setprecision(2);


  for (int k=0;k<m;k++){
  iw=1;

  while (seq_served_nodes(k,iw)!=0 && iw<2*n+2){
iw=iw+1;
//cout<<"Seq_served="<<seq_served_nodes(k,iw)<<endl;
//cout<<"k="<<k<<"\t iw="<<iw<<endl;
}


DARP_results_GLPK << "Vehicle "<<k+1<<"\n\n";
DARP_results_GLPK  << "Served nodes:\n"<<seq_served_nodes( arma::span(k,k),arma::span(0,iw-1) )<<"\n";
DARP_results_GLPK << "Service Starting time at nodes:\n"<<rel_service_start_mat( arma::span(k,k),arma::span(0,iw-1) )<<"\n";
DARP_results_GLPK << "Loads of Vehicle after leaving nodes:\n"<<rel_load_mat( arma::span(k,k),arma::span(0,iw-1) )<<"\n";
DARP_results_GLPK<<  "Ride time of the served nodes\n"<<rel_ride_time_mat.row(k)<<"\n \n \n \n";


  }


  DARP_results_GLPK.close();

//cout<<rel_ride_time_mat<<endl;



 }

    return 0;
}
