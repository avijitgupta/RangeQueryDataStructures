#include <iostream>
#include <cstdlib>
using namespace std;
#define MAXD 401 //input can vary from 0-399
int pointSet[MAXD][MAXD];
void addPoint(int x, int y);
void addYDim(int x, int y);
int searchTree(int x, int y);
int searchY(int x, int y);

void addPoint(int x, int y)
{
	while(x<MAXD)
	{
		addYDim(x, y);
		x+= (x & -x);
	}
}

void addYDim(int x, int y)
{
	while(y<MAXD)
	{
		pointSet[x][y] +=1;
		y+= (y & -y);
	}
}

int searchTree(int x, int y)
{
	int count = 0;
	while(x!= 0)
	{
		count+= searchY(x, y);
		x = x & (x-1);
	}
	return count;
}

int searchY(int x, int y)
{
	int count =0;
	while(y!=0)
	{
		count+= pointSet[x][y];
		y = y & (y-1);
	}
	return count;
}
int main()
{
	int N, R, i; 
	cin>>N>>R;

	
	pair<int,int> r_x[R];
	pair<int,int> r_y[R];
	int a, b, c, d, x , y;
	//Input pointset
	for(i=0;i<N;i++)
	{
		cin>>a>>b;
		addPoint((a+1),(b+1)); // We maintain one row to left and top as 0
	
	}
	
	for(i=0;i<R;i++)
	{
		cin>> a >> b >> c >> d;
		
		//Ensuring x1 < x2 & y1 < y2 
		// Handle degenrate case of x1 = x2 / y1 = y2?
		
		if(a<b)
		{
			r_x[i].first = a;
			r_x[i].second = b;
		} 
		else
		{
			r_x[i].first = b;
			r_x[i].second = a;
		}
		
		if(c<d)
		{
			r_y[i].first = c;
			r_y[i].second = d;
		} 
		else
		{
			r_y[i].first = d;
			r_y[i].second = c;
		}
		
	}
	
	for(i=0;i<R;i++)
	{
		int x1, x2, y1, y2;
		x1 = r_x[i].first + 1;
		x2 = r_x[i].second + 1;
		y1 = r_y[i].first + 1;
		y2 = r_y[i].second + 1;
		int bottom_right = searchTree(x2, y2);
		int bottom_left = searchTree(x2, y1-1);
		int top_left = searchTree(x1 - 1, y1 - 1);
		int top_right = searchTree(x1-1, y2);
		
		int num_points = bottom_right - top_right - bottom_left + top_left;
		cout<<num_points<<endl;
	}
}
