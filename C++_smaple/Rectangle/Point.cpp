/*
 * Point.cpp
 *
 *  Created on: 31 Mar 2020
 *      Author: kaveh
 */

#include "Point.h"
#include<iostream>
#include<list>
#include<map>
using namespace std;

Point::Point() {
	// TODO Auto-generated constructor stub
	x=0.0f;
	y=0.0f;
}
Point::Point(float xIn, float yIn) : x{xIn},y{yIn}{};
float Point::getx() const
{
	return x;
}

float Point::gety() const
{
	return y;
}
void Point::setPoint(float xIn, float yIn)
{
	x=xIn;
	y=yIn;
}
void Point::swap(double &n1,double &n2, double &n3)
{
	cout<<"Before Swapping"<<endl;
	cout<<"n1="<<n1<<"\t n2="<<n2<<"\t n3="<<n3<<endl;

	float temp;
	temp =n1;
	n1=n2;
	n2=n3;
	n3=temp;
	cout<<"After Swapping"<<endl;
	cout<<"n1="<<n1<<"\t n2="<<n2<<"\t n3="<<n3<<endl;

}
void Point::cathfunction()
{
	Point::innercath();
}

void Point::innercath()
{
	//try different exceptions
	bool catchint=false;
	bool catchchar=false;
	bool catchString=true;
	if (catchint)
	{
		throw(8);
	}

	if (catchchar)
	{
		throw("Character Err.");
	}

	if (catchString)
	{
		throw(string("String Error!"));
	}

}
bool compare(int a, int b)
	{
	    return ((int)a == (int)b);
	}
void Point::playList()
{
	list<int> lst;
	lst.push_back(1);
	lst.push_back(1);
	lst.push_back(1);
	lst.push_back(3);
	lst.push_front(0);

	for(auto it=lst.begin();it!=lst.end();it++)
	{
		if(*it==3)
		{
			lst.insert(it,4);
		}

	}
lst.unique(compare);
	for(auto it=lst.begin();it!=lst.end();it++)
	{
		cout<<*it<<" "<<endl;
	}
}
void Point::PlayMap(){
std::map<string,string> kmap;
kmap["Kaveh"]="Azizian";
kmap["Kian"]="Azizian";
kmap["Maryam"]="Azizian";
kmap["Ozra"]="Azizian";
kmap.insert(make_pair("Hosein","Azizian"));
for(auto it=kmap.begin();it!=kmap.end();it++)
{
cout<<it->first<<":"<<it->second<<endl;
}
if (kmap.find("Kian")!=kmap.end())
		{cout<<"Tapildi"<<endl;}
else
{cout<<"Tapilmadi"<<endl;}
}
Point::~Point() {
	// TODO Auto-generated destructor stub
}

