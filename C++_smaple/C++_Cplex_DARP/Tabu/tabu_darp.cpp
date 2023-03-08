#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

void getCost(double*,double*,double*,double*,int*,int [],int [],double [],double [],int [],int [],double [],double [],int,int,int,int [],int []);
double getPCost(double[],double[],int);
void getInsertion(int*,int*,double*,double*, double*, int*,int [],int [],double [],double [],int [],int [],double [],double [],int, int, int, int [], int [], int [], double []);
void removeVertex(int *, int *, int, int);
int maxOf(int [],int,int);
void copyArr(int [], int *, int);
void showResult(double *resultPtr, double *depotPtr, int p[], int n, int m, int route[], int sequence[], double, double);
void saveResult(double *resultPtr, double *depotPtr, int p[], int n, int m, int route[], int sequence[], double best_cost, double program_time, string filename);
void saveResult2(double *resultPtr, double *depotPtr, int p[], int n, int m, int route[], int sequence[], double best_cost, double program_time, string filename);

int main()
{
	/*
	-------------------- PART 1 --------------------
	    INITIALIZE

	    - get instance data from text file.
	    - calculate distances between points.
	------------------------------------------------
	*/

	string filename = "pr02";
	string filename2 = "darp/" + filename + ".txt";
	ifstream myFile (filename2.c_str());

	if (!myFile)
	{
		cerr << "File could not be opened" << endl;
	}

	// 1.1: get {m n Tk Qk L}

	int m, n, Tk, Qk, L;

	myFile >> m >> n >> Tk >> Qk >> L;

	int T[m], Q[m];

	for (int i=0; i<m; i++)
	{
		T[i] = Tk;
		Q[i] = Qk;
	}

	int N = n/2;

	// 1.2: get {x y d q e l}

	n++;

	int temp;
	double x[n], y[n], e[n], l[n];
	int d[n], q[n];

	for (int i=0; i<n; i++)
	{
		myFile >> temp >> x[i] >> y[i] >> d[i] >> q[i] >> e[i] >> l[i];
	}

	myFile.close();

	// 1.3: Calculate distances

	double distance[n][n] = {};

	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			distance[i][j] = sqrt( pow(x[i]-x[j],2)+pow(y[i]-y[j],2) );


	// 1.4: Time Window Adjustment

	{
		int N = (n-1)/2;

		for (int i=1; i<=N; i++)
		{
			e[i] = max(e[i], e[i+N]-L-d[i]);
			l[i] = min(l[i], l[i+N]-d[i]);
		}

		for (int i=N+1; i<n; i++)
		{
			e[i] = max(e[i], e[i-N]+d[i-N]);
			l[i] = min(l[i], l[i-N]+L+d[i-N]);
		}
	}

	// 1.5: Initialize Table of Result

	int best_route[n] = {};
	int best_sequence[n] = {};
	int better_route[n] = {};
	int better_sequence[n] = {};
	int current_route[n] = {};
	int current_sequence[n] = {};
	int previous_route[n] = {};
	int previous_sequence[n] = {};
	double best_cost, better_cost, current_cost, previous_cost;

	double result[n][5] = {};

	/*
	result[][0]		=>	Arrival time, A
	result[][1]		=>	Beginning of service, B
	result[][2]		=>	Departure time, D
	result[][3]		=>	Waiting time, W
	result[][4]		=>	Ride time, L
	*/

	double depot[m][2] = {};

	/*
	depot[][0]		=>	Start time (Leaving Depot)
	depot[][1]		=>	End time (Back to Depot)
	*/

	int p[n] = {};

	/*
	p[n]			=>	Number of passenger
	*/

	const int theta_max = round(7.5*log10((n-1)/2));
	const double delta_max = 0.5;
	const double lamda_max = 0.015;


	//double cost_coefficient[5] = {lamda*sqrt(N*m),1,1,1,1};

	/*
	cost_coefficient[0]	=>	coefficient for cost penalty
	cost_coefficient[0]	=>	coefficient for load violation
	cost_coefficient[1]	=>	coefficient for duration violation
	cost_coefficient[2]	=>	coefficient for time window violation
	cost_coefficient[3]	=>	coefficient for time constraint violation
	*/

	int attribute_repetition[N][m] = {};

	double cost[5] = {};

	/*
	cost[0]		=> total_cost;
	cost[1]		=> total_load_violation;
	cost[2]		=> total_time_window_violation;
	cost[3]		=> total_time_constraint_violation;
	cost[4]		=> total_duration_violation;
	*/

	/*
	-------------------- PART 2 --------------------
	    RANDOM START

	    - get instance data from text file.
	    - calculate distances between points.
	------------------------------------------------
	*/

	int iteration = 0;

	label0:
	clock_t start = clock();

	// Random Allocation

	do {
		int j=1, k=0;
		int N = (n-1)/2;

		for (int i=1; i<N+1; i++)
		{
			best_route[i] = k;
			best_route[i+N] = k;
			best_sequence[i] = j;
			best_sequence[i+N] = j+1;
			k++;

			if (k==m)
			{
				k=0;
				j+=2;
			}
		}
	} while(0);

	best_cost = 0;

	/*
	-------------------- PART 3 --------------------
	    REAL PROGRAM

	    - ??
	------------------------------------------------
	*/

	// Program Starts


	label1:

	cout << endl << best_cost;
	//try_again = 3;

	/*
	cout << endl << endl;

	for (int i=0; i<n; i++)
	{
		cout << i << "\t" << best_route[i] << "\t" << best_sequence[i] << endl;
	}

	cout << endl << best_cost << endl;
	*/


	{
		int timenow = time(NULL);
		int theta = theta_max;
		double delta = delta_max;
		double lamda = lamda_max;
		int attribute_table[n][m] = {};
		int tabo_list[theta][2] = {};
		int tabo_index = 0;
		double cost_coefficient[5] = {lamda*sqrt((n-1.0)*m/2.0),1.0,1.0,1.0,1.0};

		copyArr(best_route,&previous_route[0],n);
		copyArr(best_sequence,&previous_sequence[0],n);

		for (int iteration=0; iteration<200; iteration++)
		{

			if (iteration%20 == 19)
			{
				{
					srand(timenow=timenow+1);
					delta = rand()%21*delta_max/20;

					srand(timenow=timenow+1);
					lamda = rand()%21*lamda_max/20;
					cost_coefficient[0] = lamda*sqrt((n-1.0)*m/2.0);

					srand(timenow=timenow+1);
					theta = rand()%(theta_max+1);
				}
				/*
				for (int i=1; i<((n+1)/2); i++)
				{
					int attribute[3]={i,-1,previous_route[i]};
					removeVertex(&previous_route[0], &previous_sequence[0], n, i);
					getInsertion(&current_route[0],&current_sequence[0],&result[0][0],&depot[0][0],&distance[0][0],&p[0],previous_route,previous_sequence,x,y,d,q,e,l,m,n,L,T,Q,attribute,cost_coefficient);
					copyArr(current_route,&previous_route[0],n);
					copyArr(current_sequence,&previous_sequence[0],n);
				}*/
			}

			better_cost = 0;
			bool violation[4] = {};
			int tabo[2];

			for (int i=1; i<((n+1)/2); i++)
				for (int k=0; k<m; k++)
					if (previous_route[i] != k)
					{
						int attribute[3]={i,previous_route[i],k};
						getInsertion(&current_route[0],&current_sequence[0],&result[0][0],&depot[0][0],&distance[0][0],&p[0],previous_route,previous_sequence,x,y,d,q,e,l,m,n,L,T,Q,attribute,cost_coefficient);
						getCost(&cost[0],&result[0][0],&depot[0][0],&distance[0][0],&p[0],current_route,current_sequence,x,y,d,q,e,l,m,n,L,T,Q);
						current_cost = getPCost(cost,cost_coefficient,attribute_table[i][k]);

						//cout << endl << i << "\t" << k << "\t" << current_cost << "\t";


						if ((cost[1]+cost[2]+cost[3]+cost[4]==0) && ((cost[0] < best_cost)||(best_cost == 0)))
						{
							//cout << "1";
							copyArr(current_route,&best_route[0],n);
							copyArr(current_sequence,&best_sequence[0],n);
							best_cost = cost[0];
							goto label1;
						}
						else if (better_cost == 0)
						{
							//cout << "2";
							copyArr(current_route,&better_route[0],n);
							copyArr(current_sequence,&better_sequence[0],n);
							better_cost = current_cost;
							tabo[0] = i;
							tabo[1] = previous_route[i];
							violation[0] = cost[1];
							violation[1] = cost[2];
							violation[2] = cost[3];
							violation[3] = cost[4];
						}
						else if (current_cost < better_cost)
						{
							//cout << "3";
							bool tabo_state = 1;
							for (int j=0; j<theta; j++)
							{
								int k = tabo_index - j;

								if (k<0)
									k = k+theta_max;

								if ((tabo_list[k][0] == i) && (tabo_list[k][1] == k))
								{
									tabo_state = 0;
									break;
								}
							}

							if (tabo_state)
							{
								copyArr(current_route,&better_route[0],n);
								copyArr(current_sequence,&better_sequence[0],n);
								better_cost = current_cost;
								tabo[0] = i;
								tabo[1] = previous_route[i];
								violation[0] = cost[1];
								violation[1] = cost[2];
								violation[2] = cost[3];
								violation[3] = cost[4];
							}
						}
					}

			copyArr(better_route,&previous_route[0],n);
			copyArr(better_sequence,&previous_sequence[0],n);
			previous_cost = better_cost;
			attribute_table[tabo[0]][tabo[1]] += 1;
			tabo_list[tabo_index][0] = tabo[0];
			tabo_list[tabo_index][1] = tabo[1];
			tabo_index++;

			if (tabo_index == theta_max)
				tabo_index = 0;

			for (int i=0; i<4; i++)
			{
				if (violation[i])
					cost_coefficient[i+1] = cost_coefficient[i+1] * (1+delta);
				else
					cost_coefficient[i+1] = cost_coefficient[i+1] / (1+delta);
			}

			if (isinf(better_cost))
			{
				cout << "\ninfinity reached\n";
				break;
			}
		}

	}

	clock_t end = clock();
	double program_time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;

	if (best_cost)
	{
		getCost(&cost[0],&result[0][0],&depot[0][0],&distance[0][0],&p[0],best_route,best_sequence,x,y,d,q,e,l,m,n,L,T,Q);
		showResult(&result[0][0], &depot[0][0], p, n, m, best_route, best_sequence, best_cost, program_time);
		saveResult2(&result[0][0], &depot[0][0], p, n, m, best_route, best_sequence, best_cost, program_time, filename);
	}

	iteration++;
	if (iteration < 1)
		goto label0;

	return 0;
}

void getCost(double *costPtr, double *resultPtr,double *depotPtr, double *distancePtr, int *pPtr,int route[],int sequence[],double x[],double y[],
				int d[],int q[],double e[],double l[],int m, int n, int L, int T[], int Q[])
{
	//cout << *(costPtr+1);
	//cout << *(resultPtr+5*1+1);
	//cout << *(depotPtr+2*1+1);
	//cout << *(distancePtr+n*1+1);
	//cout << *(pPtr+1);

	int length = maxOf(sequence,0,n)+2;
	int sequenceByRoute[m][length] = {};

	for (int i=1; i<n; i++)
		sequenceByRoute[route[i]][sequence[i]]=i;

	double total_cost=0, total_load_violation=0, total_time_window_violation=0, total_time_constraint_violation=0, total_duration_violation=0;

	for (int k=0; k<m; k++)
	{
		int i_old = 0, pi_old = 0, pi;
		double Di_old, Di;

		{
			int i = sequenceByRoute[k][1];
			Di_old = e[i]-*(distancePtr+i);
			if (Di_old < 0)
				Di_old = 0;

			*(depotPtr+2*k) = Di_old;
		}

		for (int j=1; j<length; j++)
		{
			int i = sequenceByRoute[k][j];
			if (i)
			{
				double Ai,Bi,Wi,Li=0;
				Ai = Di_old + *(distancePtr+n*i_old+i);
				Bi = max(Ai,e[i]);
				Di = Bi + d[i];
				Wi = Bi - Ai;

				if (i>(n-1)/2)
					Li = Bi - *(resultPtr+5*(i-(n-1)/2)+2);

				pi = pi_old + q[i];

				double cost, load_violation, time_window_violation, time_constraint_violation;

				cost = *(distancePtr+n*i_old+i);
				load_violation = max((pi-Q[k]),0);
				time_window_violation = max((Bi-l[i]),0.0);
				time_constraint_violation = max((Li-L),0.0);

				*(resultPtr+5*i) = Ai;
				*(resultPtr+5*i+1) = Bi;
				*(resultPtr+5*i+2) = Di;
				*(resultPtr+5*i+3) = Wi;
				*(resultPtr+5*i+4) = Li;
				*(pPtr+i) = pi;

				total_cost = total_cost + cost;
				total_load_violation = total_load_violation + load_violation;
				total_time_window_violation = total_time_window_violation + time_window_violation;
				total_time_constraint_violation = total_time_constraint_violation + time_constraint_violation;
			}
			else
			{
				double cost = *(distancePtr+i_old);
				total_cost = total_cost + cost;
				*(depotPtr+2*k+1) = Di_old + *(distancePtr+n*i_old+i);
				double duration_violation = max(*(depotPtr+2*k+1)-*(depotPtr+2*k)-T[k],0.0);
				total_duration_violation = total_duration_violation + duration_violation;
				break;
			}
			i_old = i;
			Di_old = Di;
			pi_old = pi;
		}
	}

	*(costPtr) = total_cost;
	*(costPtr+1) = total_load_violation;
	*(costPtr+2) = total_time_window_violation;
	*(costPtr+3) = total_time_constraint_violation;
	*(costPtr+4) = total_duration_violation;

}

double getPCost(double cost[],double cost_coefficient[], int attribute_repetition)
{
	return( (1+ cost_coefficient[1]*attribute_repetition)*cost[0] + cost_coefficient[1]*cost[1] + cost_coefficient[2]*cost[2] + cost_coefficient[3]*cost[3] + cost_coefficient[4]*cost[4] );
}

void getInsertion(int *routePtr,int *sequencePtr,double *resultPtr,double *depotPtr, double *distancePtr, int *pPtr,int route[],int sequence[],double x[],
				double y[],int d[],int q[],double e[],double l[],int m,int n,int L,int T[],int Q[],int attribute[],double cost_coefficient[])
{
	//cout << *(resultPtr+5*1+1);
	//cout << *(depotPtr+2*1+1);
	//cout << *(distancePtr+n*1+1);
	//cout << *(pPtr+1);

	int length = 0;

	for (int i=1; i <= (n-1)/2; i++)
		if (route[i] == attribute[2])
			length++;

	length = length*2 + 4;

	int sequenceByRoute1[length] = {};
	int sequenceByRoute2[length] = {};

	for (int i=1; i<n; i++)
		if (route[i] == attribute[2])
		sequenceByRoute1[sequence[i]]=i;

	double cost_best;
	int p1_best, p2_best;
	bool state = false;

	for (int p1=1; p1<length-2; p1++)
		for (int p2=p1+1; p2<length-1; p2++)
		{
			{
				int j=0;

				for (int i=0; i<length; i++)
				{
					if (i+j == p1)
					{
						sequenceByRoute2[i+j] = attribute[0];
						j++;
					}
					if (i+j == p2)
					{
						sequenceByRoute2[i+j] = attribute[0] + (n-1)/2;
						j++;
					}
					sequenceByRoute2[i+j] = sequenceByRoute1[i];
				}
			}

			double total_cost=0, total_load_violation=0, total_time_window_violation=0, total_time_constraint_violation=0, total_duration_violation=0;

			{
				int k = attribute[2];
				int i_old = 0, pi_old = 0, pi;
				double Di_old, Di;

				{
					int i = sequenceByRoute2[1];
					Di_old = e[i]-*(distancePtr+i);
					if (Di_old < 0)
						Di_old = 0;

					*(depotPtr+2*k) = Di_old;
				}

				for (int j=1; j<length; j++)
				{
					int i = sequenceByRoute2[j];
					if (i)
					{
						double Ai,Bi,Wi,Li=0;
						Ai = Di_old + *(distancePtr+n*i_old+i);
						Bi = max(Ai,e[i]);
						Di = Bi + d[i];
						Wi = Bi - Ai;

						if (i>(n-1)/2)
							Li = Bi - *(resultPtr+5*(i-(n-1)/2)+2);

						pi = pi_old + q[i];

						double cost, load_violation, time_window_violation, time_constraint_violation;

						cost = *(distancePtr+n*i_old+i);
						load_violation = max((pi-Q[k]),0);
						time_window_violation = max((Bi-l[i]),0.0);
						time_constraint_violation = max((Li-L),0.0);

						*(resultPtr+5*i) = Ai;
						*(resultPtr+5*i+1) = Bi;
						*(resultPtr+5*i+2) = Di;
						*(resultPtr+5*i+3) = Wi;
						*(resultPtr+5*i+4) = Li;
						*(pPtr+i) = pi;

						total_cost = total_cost + cost;
						total_load_violation = total_load_violation + load_violation;
						total_time_window_violation = total_time_window_violation + time_window_violation;
						total_time_constraint_violation = total_time_constraint_violation + time_constraint_violation;
					}
					else
					{
						double cost = *(distancePtr+i_old);
						total_cost = total_cost + cost;
						*(depotPtr+2*k+1) = Di_old + *(distancePtr+n*i_old+i);
						double duration_violation = max(*(depotPtr+2*k+1)-*(depotPtr+2*k)-T[k],0.0);
						total_duration_violation = total_duration_violation + duration_violation;
						break;
					}
					i_old = i;
					Di_old = Di;
					pi_old = pi;
				}
			}

			double cost[5] = {total_cost, total_load_violation, total_time_window_violation, total_time_constraint_violation, total_duration_violation};
			double cost_now = getPCost(cost,cost_coefficient,0);


			if (state == false)
			{
				p1_best = p1;
				p2_best = p2;
				cost_best = cost_now;
				state = true;
			}
			else if ( cost_now < cost_best)
			{
				p1_best = p1;
				p2_best = p2;
				cost_best = cost_now;
			}
		}

	int initial_p1 = sequence[attribute[0]];
	int initial_p2 = sequence[attribute[0]+(n-1)/2];

	for (int i=1; i<n; i++)
		if (route[i] == attribute[1])
		{
			if (sequence[i] > initial_p2)
			{
				*(routePtr+i) = route[i];
				*(sequencePtr+i) = sequence[i]-2;
			}
			else if (sequence[i] > initial_p1)
			{
				*(routePtr+i) = route[i];
				*(sequencePtr+i) = sequence[i]-1;
			}
			else
			{
				*(routePtr+i) = route[i];
				*(sequencePtr+i) = sequence[i];
			}
		}
		else if (route[i] == attribute[2])
		{
			if (sequence[i] > p2_best-2)
			{
				*(routePtr+i) = route[i];
				*(sequencePtr+i) = sequence[i]+2;
			}
			else if (sequence[i] > p1_best-1)
			{
				*(routePtr+i) = route[i];
				*(sequencePtr+i) = sequence[i]+1;
			}
			else
			{
				*(routePtr+i) = route[i];
				*(sequencePtr+i) = sequence[i];
			}
		}
		else
		{
			*(routePtr+i) = route[i];
			*(sequencePtr+i) = sequence[i];
		}

	int p1 = attribute[0];
	int p2 = attribute[0] + (n-1)/2;

	*(routePtr+p1) = attribute[2];
	*(routePtr+p2) = attribute[2];
	*(sequencePtr+p1) = p1_best;
	*(sequencePtr+p2) = p2_best;
}

void removeVertex(int *routePtr, int *sequencePtr, int n, int p1)
{
	int p2 = p1 + (n-1)/2;
	int p_route = *(routePtr+p1);
	int p1_sequence = *(sequencePtr+p1);
	int p2_sequence = *(sequencePtr+p2);

	for (int i=1; i<n; i++)
		if (*(routePtr+i) == p_route)
		{
			if (*(sequencePtr+i) > p2_sequence)
			{
				*(sequencePtr+i) = *(sequencePtr+i)-2;
			}
			else if (*(sequencePtr+i) > p1_sequence)
			{
				*(sequencePtr+i) = *(sequencePtr+i)-1;
			}
		}

	*(routePtr+p1) = -1;
	*(routePtr+p2) = -1;
	*(sequencePtr+p1) = 1;
	*(sequencePtr+p2) = 2;
}

int maxOf(int arr[],int a, int b)
{
	int maxnum = arr[a];
	for (int i=a; i<b; i++)
		if (arr[i]>maxnum)
			maxnum = arr[i];
	return(maxnum);
}

void copyArr(int arr1[], int *arr2, int n)
{
	for (int i=0; i<n; i++)
		*(arr2+i) = arr1[i];
}

void showResult(double *resultPtr, double *depotPtr, int p[], int n, int m, int route[], int sequence[], double best_cost, double program_time)
{
	int reqPerRoute[m] = {};

	for (int i=1; i<n; i++)
		for (int j=0; j<m; j++)
			if (route[i] == j)
			{
				reqPerRoute[j] = reqPerRoute[j] + 1;
				break;
			}

	for (int j=0; j<m; j++)
	{
		cout << endl << "\t---------- ROUTE " << j << " ----------" << endl;
		cout << endl << "\t\tA\tB\tD\tW\tL\tp";
		cout << endl << "\t0\t" << *(depotPtr + 2*j) << endl;

		for (int k=1; k<reqPerRoute[j]+1; k++)
		{
			for (int i=1; i<n; i++)
				if ((route[i] == j) && (sequence[i] == k))
				{
					cout << "\t" << i;

					for (int l=0; l<5; l++)
						cout << "\t" << *(resultPtr + 5*i + l);

					cout << "\t" << p[i] << endl;
					break;
				}
		}

		cout << "\t0\t" << *(depotPtr + 2*j + 1) << endl;
	}
	cout << endl << endl << "The minimum cost is " << best_cost;
	cout << endl << "The program duration is " << program_time << " ms";
}

void saveResult(double *resultPtr, double *depotPtr, int p[], int n, int m, int route[], int sequence[], double best_cost, double program_time, string filename)
{
	double previous_cost = 0;
	string filename2 = "benchmark/" + filename;

	ifstream myFile (filename2.c_str());

	if (myFile)
	{
		myFile >> previous_cost;
	}

	myFile.close();

	if ((best_cost < previous_cost)||(previous_cost==0))
	{
		if (best_cost < previous_cost)
			remove(filename2.c_str());

		ofstream myFile(filename2.c_str());

		myFile << best_cost <<  endl;

		int reqPerRoute[m];

		for (int i=1; i<n; i++)
			for (int j=0; j<m; j++)
				if (route[i] == j)
				{
					reqPerRoute[j] = reqPerRoute[j] + 1;
					break;
				}

		for (int j=0; j<m; j++)
		{
			myFile << endl << "\t---------- ROUTE " << j << " ----------" << endl;
			myFile << endl << "\t\tA\tB\tD\tW\tL\tp";
			myFile << endl << "\t0\t" << *(depotPtr + 2*j) << endl;

			for (int k=1; k<reqPerRoute[j]+1; k++)
			{
				for (int i=1; i<n; i++)
					if ((route[i] == j) && (sequence[i] == k))
					{
						myFile << "\t" << i;

						for (int l=0; l<5; l++)
							myFile << "\t" << *(resultPtr + 5*i + l);

						myFile << "\t" << p[i] << endl;
						break;
					}
			}

			myFile << "\t0\t" << *(depotPtr + 2*j + 1) << endl;
		}
		myFile << endl << endl << "The minimum cost is " << best_cost;
		myFile << endl << "The program duration is " << program_time << " ms";
		myFile.close();
	}
}

void saveResult2(double *resultPtr, double *depotPtr, int p[], int n, int m, int route[], int sequence[], double best_cost, double program_time, string filename)
{
	for (int i=1; ; i++)
	{
		char temp[5];
		itoa(i,temp,10);
		string filename2 = "log/" + filename + "_" + temp + ".txt";

		ifstream myFile (filename2.c_str());

		if (myFile)
		{
			myFile.close();
		}
		else
		{
			myFile.close();
			{
				ofstream myFile(filename2.c_str());

				myFile << "Solution found!" <<  endl << endl;
				myFile << "Total cost is " <<  best_cost << endl << endl << endl;

				int reqPerRoute[m];

				for (int i=1; i<n; i++)
					for (int j=0; j<m; j++)
						if (route[i] == j)
						{
							reqPerRoute[j] = reqPerRoute[j] + 1;
							break;
						}

				for (int j=0; j<m; j++)
				{
					myFile << "Vehicle " <<  j+1 << endl << endl;

					// 1) Served nodes
					myFile << "Served nodes: " << endl;
					myFile << "     0";

					for (int k=1; k<reqPerRoute[j]+1; k++)
					{
						for (int i=1; i<n; i++)
							if ((route[i] == j) && (sequence[i] == k))
							{
								myFile << "     " <<i;
								break;
							}
					}

					myFile << "     0" << endl << endl;

					// 2) Service starting time at nodes
					myFile << "Service starting time at nodes: " << endl;
					myFile << "     " << *(depotPtr + 2*j);

					for (int k=1; k<reqPerRoute[j]+1; k++)
					{
						for (int i=1; i<n; i++)
							if ((route[i] == j) && (sequence[i] == k))
							{
								myFile << "     " <<*(resultPtr + 5*i + 1);
								break;
							}
					}

					myFile << "     " << *(depotPtr + 2*j + 1) << endl << endl;

					// 3) Loads of vehicles after leaving nodes
					myFile << "Loads of vehicles after leaving nodes: " << endl;
					myFile << "     0";

					for (int k=1; k<reqPerRoute[j]+1; k++)
					{
						for (int i=1; i<n; i++)
							if ((route[i] == j) && (sequence[i] == k))
							{
								myFile << "     " <<p[i];
								break;
							}
					}

					myFile << "     0" << endl << endl;

					// 4) Ride time of the served nodes
					myFile << "Ride time of the served nodes: " << endl;
					myFile << "     0";

					for (int k=1; k<reqPerRoute[j]+1; k++)
					{
						for (int i=1; i<n; i++)
							if ((route[i] == j) && (sequence[i] == k))
							{
								myFile << "     " <<*(resultPtr + 5*i + 4);
								break;
							}
					}

					myFile << "     0" << endl << endl;


					// 5) Waiting time of the served nodes
					myFile << "Waiting time of the served nodes: " << endl;
					myFile << "     " << *(depotPtr + 2*j);

					for (int k=1; k<reqPerRoute[j]+1; k++)
					{
						for (int i=1; i<n; i++)
							if ((route[i] == j) && (sequence[i] == k))
							{
								myFile << "     " <<*(resultPtr + 5*i + 3);
								break;
							}
					}

					myFile << "     0" << endl << endl << endl << endl;


				}
				myFile.close();
			}
			break;
		}
	}
}
