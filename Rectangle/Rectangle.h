/*
 * Rectangle.h
 *
 *  Created on: 31 Mar 2020
 *      Author: kaveh
 */

#ifndef RECTANGLE_H_
#define RECTANGLE_H_
#include"Point.h"
#include<iostream>
#include <stdexcept>
using namespace std;
class Rectangle {
public:
	Rectangle();
	Rectangle(Point BRIn,Point BLIn,Point TRIn,Point TLIn) ;
	float width() const;
	float Length() const;
	float Area()const;
	float Circum() const;
	virtual float diag() const;
	virtual ~Rectangle();
protected:
	Point BL=Point(0.0f,0.0f);
	Point BR=Point(1.0f,0.0f);
	Point TL=Point(0.0f,1.0f);
	Point TR=Point(1.0f,1.0f);


};

#endif /* RECTANGLE_H_ */
