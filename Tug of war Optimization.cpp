#include <iostream>
#include <stdlib.h>
#include <stdio.h> 
#include <time.h>
#include <cmath>
using namespace std;
/* run this program using the console pauser or add your own getch, system("pause") or input loop */

// generate uniformly distributed random number
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double fitness(double x, double y){
	double a = (4-2.1*(pow(x,2))+(pow(x,4))/3)*(pow(x,2))+x*y+(-4 + 4*(pow(y,2)))*(pow(y,2));
    return a;
}

void swap(double *xp, double *yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}
 //BubbleSort is the simplest sorting function that
 // works by repeatedly swapping the fitness elements and 
 // crossponding postions of the teams  if they are in wrong order
 
void bubbleSort(double fit[], double xarr[], double yarr[], int n)
{
	int i, j;
	for (i = 0; i < n-1; i++){
		for (j = 0; j < n-i-1; j++){
			if (fit[j] > fit[j+1]){
				swap(&fit[j], &fit[j+1]);
				swap(&xarr[j], &xarr[j+1]);
				swap(&yarr[j], &yarr[j+1]);
			}
		}
	}              
}

double max(double arr[], int n){
	double max = arr[0];
	for(int i=1;i<n;i++)
		if(arr[i]>max)
			max = arr[i];
	return max;
}

double min(double arr[], int n){
	double min = arr[0];
	for(int i=1;i<n;i++)
		if(arr[i]<min)
			min = arr[i];
	return min;
}


//create a one dimensional array and use a random number generator to populate the array

double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do{
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
	 n2_cached = 0;
     return n2*stddev + mean;
    }
}

int main(int argc, char** argv) {
	
	srand(time(NULL));
	double minxrange=-2.0;
	double maxxrange=2.0;
	double minyrange=-3.0;
	double maxyrange=3.0;
	
	const int N=15;//N random candidate solution, population of the teams
	int n=2; //(x,y) of fitness function
	
	double ftns[N]; //fitness
	double weight[N]; //weight of the teams
	
	double teamsx[N];
	double teamsy[N];

	for (int i=0; i<N; i++){
	    teamsx[i] = fRand(minxrange,maxxrange);
	    teamsy[i] = fRand(minyrange,maxyrange);
	}
	
	double Ms = 1;//the coefficient of static (taken as unity Us)
	double Mk = 0.9; // 0.1 => 1 (coefficient of kinematic friction Uk)
	double dt = 5; //Delta time
	double alpha = 0.99; // [0.9 , 0.99]
	double beta = 0.05; // (0,1]

	int k=1;
	while(k <= 1000){
		for (int i=0; i<N; i++){
		    ftns[i] = fitness(teamsx[i],teamsy[i]);
		}
		
		bubbleSort(ftns,teamsx,teamsy,N);//call the function BubbleSort arranges the league of teams 
		for (int i=0; i<N; i++){
			double diffrence = (ftns[0]-ftns[N-1]);//ftns[0]=fitnessBest & ftns[N-1]=fitnessWorst 
			if(diffrence != 0)
				weight[i] = (ftns[i]-ftns[N-1])/diffrence + 1; //ftns[i] fitness value for i particle
		}
		for (int i=0; i<N; i++){
			double dxi = 0; // 
			double dyi = 0;
			for (int j=0; j<N; j++){
	            if (weight[j]*Ms > weight[i]*Ms){
	            	double frij = weight[j]*Ms - weight[i]*Mk; //pull force pour j - resistance pour i
	            	//the gravitational acceleration constant
	                double gxij = teamsx[j] - teamsx[i]; //teamsx[j]&teamsy[j] are the position vectors for candidate solution
	                double gyij = teamsy[j] - teamsy[i];
	                //axji & ayij is the cordinations of acceleration of the team
	                double axij = (frij/(weight[i]*Mk))*gxij;//frij is the pulling force between team i towards team j
	                double ayij = (frij/(weight[i]*Mk))*gyij;//weight[i]*Mk friction force 
	                //the random portion of the search space traveled by team
	                //alpha is to gradually decrease the random portionof the team's movement 
	                //beta is a scaling factor which can be chosen from the interval (0,1]
	                double dxij = 0.5*axij*pow(dt,2) + alpha*beta*(max(teamsx,N)-min(teamsx,N))*rand_normal(1,n);
	                double dyij = 0.5*ayij*pow(dt,2) + alpha*beta*(max(teamsy,N)-min(teamsy,N))*rand_normal(1,n);
	               //the new position of the team i at the end of the end of the k iteartion
				    
					dxi=dxi+dxij;
	                dyi=dyi+dyij;
				}
			}
			////the sum of the spaces that the team moved 
			teamsx[i]=teamsx[i]+dxi;
			//if the candidate solution leave the search space 
        	if(teamsx[i]<minxrange || teamsx[i]>maxxrange)
	            teamsx[i] = teamsx[0] + (rand_normal(0,1)/k)*(teamsx[0] - teamsx[i]);
	        teamsy[i]=teamsy[i]+dyi;
	        if(teamsy[i]<minyrange || teamsy[i]>maxyrange)
	            teamsy[i] = teamsy[0] + (rand_normal(0,1)/k)*(teamsy[0] - teamsy[i]);
	            //the side constraint handling technique to regenerate violating variables
	       	double newteamFit = fitness(teamsx[i],teamsy[i]);
			if(newteamFit<ftns[N-1])
            	ftns[N-1] = newteamFit;
		}
		alpha = alpha * 0.99;
		cout << min(ftns,N) << endl;
		k++;
	}
//	for (int i=0; i<N; i++){
//		cout << teamsx[i] << "\t\t" << teamsy[i] << "\t\t" << fit[i] << "\t\t" << w[i] << endl;
//	}
		return 0;
}
