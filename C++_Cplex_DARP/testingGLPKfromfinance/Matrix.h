//
// Matrix.h
 
#ifndef __FinancialSamples__Matrix__
#define __FinancialSamples__Matrix__
 
#include <vector>
 
class Matrix {
public:
typedef std::vector<double> Row;
 
Matrix(int size, int size2);
Matrix(int size);
Matrix(const Matrix &s);
~Matrix();
Matrix &operator=(const Matrix &s);
 
void transpose();
double trace();
void add(const Matrix &s);
void subtract(const Matrix &s);
void multiply(const Matrix &s);
 
Row & operator[](int pos);
private:
std::vector<Row> m_rows;
};
// free operators
//
Matrix operator+(const Matrix &s1, const Matrix &s2);
Matrix operator-(const Matrix &s1, const Matrix &s2);
Matrix operator*(const Matrix &s1, const Matrix &s2);
 
#endif /* defined(__FinancialSamples__Matrix__)
 

