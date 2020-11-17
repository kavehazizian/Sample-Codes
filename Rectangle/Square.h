/*
 * Square.h
 *
 *  Created on: 1 Apr 2020
 *      Author: kaveh
 */

#ifndef SQUARE_H_
#define SQUARE_H_
#include"Rectangle.h"
class Square :public Rectangle{
public:
	Square();
	Square(Point p1,Point p2,Point p3,Point p4 ):Rectangle(p1,p2,p3,p4){};
	float side() const;
	Point setBLcorner(float a, float b);
	Point setBRcorner(float a, float b);
	void settopcorners(Point BL, Point BR);
	float diag() const;
	virtual ~Square();
private:

};

#endif /* SQUARE_H_ */
