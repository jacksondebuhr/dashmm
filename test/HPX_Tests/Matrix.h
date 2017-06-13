#ifndef Matrix_h_
#define Matrix_h_

#include <vector>
#include <iostream>
using namespace std;


class Matrix{
private:
	//REPRESENTATION
	int n;//rows
	int m;//cols
	vector<vector<double> > A; 
public:
	Matrix();
	Matrix(int n_,int m_,double val = 0);
	Matrix(const Matrix& B);
	void operator=(const Matrix& B);

	int get_height()const;
	int get_width()const;
	void set_entry(int i,int j, double a);
	double get_entry(int i,int j)const;
	vector<double> get_col(int j)const;
	void set_col(int j, const vector<double>& v);
	vector<double> toVector()const;
	void print()const;

	//matrix and vector operations
	vector<double> mv_mult(const vector<double>& x)const;
	Matrix rh_multiply(const Matrix& B)const;
	Matrix operator-()const;
	Matrix operator+(const Matrix& B)const;
	Matrix operator-(const Matrix& B)const;
	double norm_2()const;
	double dot(const Matrix& v)const;
	Matrix scalar_multiply(double a)const;

	//jacobi method functions
	Matrix get_Dinv()const;
	Matrix get_LplusU()const;

	//GMRES method functions
	Matrix get_subMatrix(int p, int q)const;
	
};
#endif//Matrix


