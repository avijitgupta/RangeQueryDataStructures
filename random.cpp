#include <iostream>
#include <cstdlib>
#include <ctime>
#define POINTS 800
#define RANGES 50000
using namespace std;
int main()
{
	cout<<POINTS<<" "<<RANGES<<endl;
	srand(time(0));
	for(int i = 0;i<POINTS;i++)
	{
		cout<<rand()%100000<<" "<< rand()%100000<<endl;
	}

	for(int i = 0;i<RANGES;i++)
	{
		cout<<rand()%100000<<" "<< rand()%100000<<" "<<rand()%100000<<" "<<rand()%100000<<endl;
	}



}
