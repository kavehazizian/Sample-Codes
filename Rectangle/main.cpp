/*
 * main.cpp
 *
 *  Created on: 31 Mar 2020
 *      Author: kaveh
 */

#include <iostream>
#include "Point.h"
#include"Rectangle.h"
#include"square.h"
#include<vector>
using namespace std;

int main()
{
// lap avalda baxim goorim	
/*Point p1(2,4);
cout<<"Xcoordinate is"<<p1.getx()<<endl;
cout<<"ycoordinate is"<<p1.gety()<<endl;
cout<<"Now changing coordinate"<<endl;
p1.setPoint(5.1,6.4);
cout<<"Xcoordinate is"<<p1.getx()<<endl;
cout<<"ycoordinate is"<<p1.gety()<<endl;
  cout << "Ye lk siximi world!" << endl;*/

	Point BL(4,5);
	Point BR(6,5);
	Point TL(4,8);
	Point TR(6,8);

	Rectangle R(BR,BL,TR,TL);
	cout<<"Area is "<<R.Area()<<endl;
	cout<<"Circum is "<<R.Circum()<<endl;
//std::unique_ptr<Rectangle> p(new Rectangle());
	Point swp(0,0);
	double n1=3.1;
	double n2=4.5;
	double n3=6.7;
	swp.swap(n1,n2,n3);
Square s1;
//s1.setBLcorner(4,6);
//s1.setBRcorner(9,6);
// for testing the submodule
s1.settopcorners(Point(4,6),Point(9,6));

cout<<"Side is "<<s1.side()<<"\t"<<s1.Circum()<<"\t"<<s1.Area()<<endl;
Rectangle *ptr=new Rectangle;
cout<<ptr->diag()<<endl;

Rectangle *ptr2=new Square;
cout<<ptr2->diag()<<endl;
Square s2(Point(0,0),Point(3,0),Point(0,3),Point(3,3));
cout<<s2.Area() << "\t" << s2.diag()<<endl;

Point Err(0,0);
try{
Err.cathfunction();
}
catch(int e)
{
cout<<"Int Err."<<e<<endl;
}
catch(char const* e )
{
cout<<"Char Err."<<e<<endl;
}

catch(string& e)
{
cout<<"String Err."<<e<<endl;
}

Err.playList();
Err.PlayMap();
MiniPoint Mp(1,"kl");
cout<<Mp.val<<"\t"<<Mp.key<<endl;
vector<MiniPoint> vec;
vec.push_back(MiniPoint(2,"hj"));
vec.push_back(MiniPoint(3,"fgf"));
vec.push_back(MiniPoint(4,"er"));
vec.push_back(MiniPoint(5,"tyu"));
for(auto ip=vec.begin();ip!=vec.end();ip++)
{
cout<<*ip <<endl;
}
	return 0;
}


