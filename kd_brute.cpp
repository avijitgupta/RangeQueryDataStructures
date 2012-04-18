#include <iostream>
#include <cstdlib>
#include <algorithm>
#define EPSILON 0.00000000001
using namespace std;
int double_gt(double a, double b);
int double_equal(double a, double b);
int double_lt(double a, double b);
int double_equal(double a, double b)
{
	if( ((a-b) < EPSILON && (a-b) >= 0) || ((b-a) < EPSILON && (b-a) >= 0) )
		return 1;
	return 0;
	     
}

int double_lt(double a, double b)
{
	if((b - a) > EPSILON) return 1;
	return 0;
}


int double_gt(double a, double b)
{
	if((a - b) > EPSILON) return 1;
	return 0;
}

int main()
{
	
	int N, R, i,j; 
	cin>>N>>R;
	pair<double,double> p_sx[N], p_sy[N];
	
	pair<double,double> r_x[R];
	pair<double,double> r_y[R];
	double a, b, c, d;
	//Input pointset
	for(i=0;i<N;i++)
	{
		cin>>p_sx[i].first>>p_sx[i].second;
		p_sy[i].first = p_sx[i].first;
		p_sy[i].second = p_sx[i].second;
	
	}
	
	for(j=0;j<R;j++)
	{
		cin>> a >> b >> c >> d;
		
		//Ensuring x1 < x2 & y1 < y2 
		// Handle degenrate case of x1 = x2 / y1 = y2?
		
		if(double_lt(a,b))
		{
			r_x[j].first = a;
			r_x[j].second = b;
		} 
		else
		{
			r_x[j].first = b;
			r_x[j].second = a;
		}
		
		if(double_lt(c,d))
		{
			r_y[j].first = c;
			r_y[j].second = d;
		} 
		else
		{
			r_y[j].first = d;
			r_y[j].second = c;
		}
		int count =0;
		for(i=0;i < N; i++)
		{
			if( (double_equal(p_sx[i].first,r_x[j].first) || double_gt(p_sx[i].first, r_x[j].first)) && 
			    (double_equal(p_sx[i].first,r_x[j].second) || double_lt(p_sx[i].first, r_x[j].second)) &&
			    (double_equal(p_sx[i].second,r_y[j].first) || double_gt(p_sx[i].second, r_y[j].first)) && 
			    (double_equal(p_sx[i].second,r_y[j].second) || double_lt(p_sx[i].second, r_y[j].second)) )
				{	count++;	    
					cout<<p_sx[i].first<<" "<<p_sx[i].second<<endl;
				}
		}

		cout<<count<<endl;
		
	}

}
