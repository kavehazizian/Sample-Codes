/*
 * Square.cpp
 *
 *  Created on: 1 Apr 2020
 *      Author: kaveh
 */

#include "Square.h"
#include<math.h>
#include<iostream>
using namespace std;

Square::Square() {
	// TODO Auto-generated constructor stub
	cout<<"Square constructor"<<endl;
}
float Square::side() const
{
	float a =Rectangle::Length();
	float b =Rectangle::width();
	if(a!=b)
	{
		throw std::invalid_argument("This cannot be a square!");
	}
	return a;
}
Point Square::setBLcorner(float a, float b)
{
	this->BL.x=a;
	this->BL.y=b;
	return this->BL;
}
Point Square::setBRcorner(float a, float b)
{
	this->BR.x=a;
	this->BR.y=b;
	return this->BR;
}
void Square::settopcorners(Point BLin, Point BRin)
{
	this->BL=BLin;
	BR=BRin;
	this->TL.x=BLin.x;
	this->TR.x=BRin.x;
	float length=abs(BLin.x-BRin.x);
	this->TL.y=BLin.y+length;
	this->TR.y=TL.y;
}

float Square::diag() const
{
cout<<"Square diag here"<<endl;
return(sqrt(2)*this->side());

}
Square::~Square() {
	// TODO Auto-generated destructor stub
}

