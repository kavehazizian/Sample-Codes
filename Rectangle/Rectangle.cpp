/*
 * Rectangle.cpp
 *
 *  Created on: 31 Mar 2020
 *      Author: kaveh
 */

#include "Rectangle.h"
#include"Point.h"
#include <math.h>

Rectangle::Rectangle() {
	// TODO Auto-generated constructor stub
	//BL=Point(1,1);

}
Rectangle::Rectangle(Point BLIn,Point BRIn,Point TLIn,Point TRIn) : BL{BLIn.x,BLIn.y},BR{BRIn.x,BRIn.y},TL{TLIn.x,TLIn.y},TR{TRIn.x,TRIn.y}
{
	//Check out the correct coordinates
	//try{
	if ((BRIn.x != TRIn.x)||(BRIn.y != BLIn.y)|| (BLIn.x != TLIn.x)||(TRIn.y != TLIn.y))
		throw std::invalid_argument( "The coordinates cannot be a rectangle" );
	//}
	//catch(const std::exception& e)
	//{
	//	cout<<"The coordinates cannot be a rectangle222"<<endl;
	//}

};

float Rectangle::Length()const
{
	float l=sqrt(pow((BL.x-BR.x),2)+pow((BL.y-BR.y),2));
	return l;
}
float Rectangle::width() const
{
	float w=sqrt(pow((BL.x-TL.x),2)+pow((BL.y-TL.y),2));
	return w;
}
float Rectangle::Area() const
{
	return this->width()*this->Length();
}

float Rectangle::Circum() const
{
	return 2*(this->width()+this->Length());
}
float Rectangle::diag() const

{
	cout<<"Calling Rectangle diag"<<endl;
	return sqrt(pow(this->Length(),2)+pow(this->width(),2));
}
Rectangle::~Rectangle() {
	// TODO Auto-generated destructor stub
}

