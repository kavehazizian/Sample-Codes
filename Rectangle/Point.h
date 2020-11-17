/*
 * Point.h
 *
 *  Created on: 31 Mar 2020
 *      Author: kaveh
 */

#ifndef POINT_H_
#define POINT_H_
#include<iostream>
using namespace std;
class Point {
public:
	Point();
	Point(float xIn=0.0, float yIn=0.0);
	virtual ~Point();
	void setPoint(float x, float y);
	void swap(double &n1,double &n2, double &n3);
	float getx() const;
	float gety() const;
	float x;
	float y;
	void cathfunction();
	void innercath();
	void playList();
	void PlayMap();
};
class MiniPoint
{
public:
	int val;
	string key;
	MiniPoint(int val1=0, string key1=""): val(val1), key(key1){};
	MiniPoint(const MiniPoint& other)
	{
		*this=other;
	}
	void print()
	{
		cout<<"Key="<<key<<"\t val"<<val<<endl;

	}
	friend ostream &operator<<(ostream& out, const MiniPoint& MP)
	{
		out<<MP.key<<":"<<MP.val;
		return out;

	}

	};
#endif /* POINT_H_ */
