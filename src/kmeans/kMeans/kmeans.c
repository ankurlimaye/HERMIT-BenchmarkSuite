#include <iostream>
#include <stdlib.h>
//#include <conio.h>
#include <math.h>
#include <vector>

using namespace std;

int min(int arr[], int maxIndex)
{
	int min=100;
	for(int i=0;i<maxIndex;i++)
	{
		if(arr[i]<min)
		min=arr[i];
	}
	return min;
}
int indexOf(int number,int arr[], int maxIndex)
{
	int index;
	for(int i=0;i<maxIndex;i++)
	{
		if(number==arr[i])
		{
		index=i;
		break;	
		}
	}
	return index;
}
int mean(vector<int> vc )
{
int sum=0;
for(int i=0;i<vc.size();i++)
sum=sum+vc[i];
return sum/vc.size();
}
void show(vector<int> vc )
{
for(int i=0;i<vc.size();i++){
cout<<vc[i]<<",";	
}
}
bool isEqual(int arr1[], int arr2[], int maxIndex){
	for(int i=0;i<maxIndex;i++)
	{
		if(arr1[i]!=arr2[i])
		return false;
	}
	return true;
}
int main()
{
	
	int noOfItems;
	int k;	
	noOfItems = 100;
	k = 10;
	int cluster[k];
	int oldCluster[k];
	int objects[noOfItems];
	int row[k];
	vector< vector<int> > groups; 
	
	for(int i=0;i<noOfItems;i++) 
	{
		objects[i] = rand() % 50;
		if(i<k) 
		cluster[i]=objects[i];
	}
	for(int i=0;i<k;i++)
	{
	vector<int> newGroup;
	groups.push_back(newGroup);
	}
	int iter =1;
	do
	{
	for(int i=0;i<noOfItems;i++)
	{
		for(int j=0;j<k;j++){
		row[j] = abs(cluster[j]-objects[i]); 
		}		
	    groups[indexOf(min(row,k),row,k)].push_back(objects[i]); 
	}
	
	for(int j=0;j<k;j++)
	{
    	if(!groups[j].empty())
		{
		oldCluster[j]=cluster[j]; 
		cluster[j] = mean(groups[j]); 	
		}
	}
	if(!isEqual(oldCluster,cluster,k))
	{
		for(int i=0;i<k;i++)
		groups[i].clear();
	}
	iter++;	
	}while(!isEqual(oldCluster,cluster,k)); 
	cout<<"nn"; 
	for(int i=0;i<k;i++) 
	{
		cout<<"C"<<(i+1)<<" : "<<cluster[i]<<endl;
	}
	for(int i=0;i<k;i++)
	{
		cout<<"nnGroup "<<(i+1)<<" : n"<<endl;
		show(groups[i]);
	}
	cout<<"nnNumber of Iterations "<<iter<<endl;
	return 0;
}
