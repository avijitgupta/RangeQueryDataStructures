#include <iostream>
#include <cstdlib>
using namespace std;
#define LOW_X 0
#define LOW_Y 0
#define HIGH_X 40
#define HIGH_Y 40
using namespace std;
int main()
{
	int i, j, k, l;
	int count = (HIGH_X/2) * (HIGH_Y/2);
	cout<<count;
	count = 0;
	for(i = 1;i< HIGH_X ; i+=2)
	{
		for(j = i+2; j<HIGH_Y; j+=1)
		{
			for(k=1;k<HIGH_X;k+=2)
			{
				for(l=k+2;l<HIGH_Y; l+=1)
				{
					count++;
				}
			}
		}
	}
	cout<<" "<<count<<endl;
	for(i=0;i<HIGH_X; i+=2)
	{
		for(j = 0; j<HIGH_Y; j+=2)
		{
			cout<<i<<"\t"<<j<<endl;
		}
	}
	
	for(i = 1;i< HIGH_X ; i+=2)
	{
		for(j = i+2; j<HIGH_Y; j+=1)
		{
			for(k=1;k<HIGH_X;k+=2)
			{
				for(l=k+2;l<HIGH_Y; l+=1)
				{
					cout<<i<<" "<<j<<" "<<k<<" "<<l<<endl;
				}
			}
		}
	}
}
