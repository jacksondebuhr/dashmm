#include "Matrix.h"
#include <math.h> //sqrt

Matrix::Matrix(){}

Matrix::Matrix(int n_,int m_,double val){
	n = n_; m = m_;
	//vector<vector<double> > A;
	for(int i=0;i<n_;++i){
		vector<double> row;
		for(int j=0;j<m_;++j){
			row.push_back(val);
		}
		A.push_back(row);
	}
}

Matrix::Matrix(const Matrix& B){
	n = B.n; m = B.m;
	for(int i=0;i<n;++i){
		vector<double> row;
		for(int j=0;j<m;++j){
			row.push_back(B.get_entry(i,j));
		}
		A.push_back(row);
	}
}

void Matrix::operator=(const Matrix& B){
	n = B.n; 
	m = B.m;
	A = B.A;
}

int Matrix::get_height()const{
	return n;
}

int Matrix::get_width()const{
	return m;
}

void Matrix::set_entry(int i,int j, double a){
	A[i][j] = a;
}

double Matrix::get_entry(int i,int j)const{
	double d = A[i][j];
	return d;
}

vector<double> Matrix::get_col(int j)const{
	vector<double> jth_col(n);//preallocate
	for(int i=0;i<n;++i){
		jth_col[i] = A[i][j];
	}
	return jth_col;
}

void Matrix::set_col(int j, const vector<double>& v){
	//replace less than and up to the entire col
	if(v.size() <= n && j<m){
		for(int i=0;i<v.size();++i){
			A[i][j] = v[i];
		}
	}else{
		cerr<<"set_col: Dimensions do not match!"<<endl;
	}
}

vector<double> Matrix::toVector()const{
	vector<double> v;
	if(m == 1){
		for(int i=0;i<n;++i){
			v.push_back(A[i][0]);
		}
	}else if(n == 1){
		for(int j=0;j<n;++j){
			v.push_back(A[0][j]);
		}
	}else{
		cerr<<"Not a 1D Matrix"<<endl;
	}
	return v;
}

void Matrix::print()const{
	//cout<<"in print"<<endl;
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			cout<<A[i][j]<<" ";
		}
		cout<<endl;
	}
}

/*/////////////////////////////////////////////////////////////////////////
			Matrix Operations
/////////////////////////////////////////////////////////////////////////*/

//USE DASHMM TO ACHIEVE O(n) EVENTUALLY
vector<double> Matrix::mv_mult(const vector<double>& x)const{
	int q = x.size();
	if(m != q){
		cerr<<"mv_mult: Dimensions do not match!"<<endl;
	}
	vector<double> y;
	for(int i=0;i<n;++i){
		double e = 0;
		for(int j=0;j<m;++j){
			e = e + A[i][j]*x[j];
		}	
		y.push_back(e);
	}
	return y;
}

//computes (THIS)*B
Matrix Matrix::rh_multiply(const Matrix& B)const{
	//THIS:nxm  B:pxq
	//cout<<"in rh_multiply"<<endl;
	int p = B.get_height();
	if(m != p){
		cerr<<"rh_multiply: Dimensions do not match!"<<endl;
	}
	int q = B.get_width();
	Matrix C(n,q);
	//cout<<n<<" x "<<m<<"  "<<p<<" x "<<q<<endl;
	for(int i=0;i<n;++i){
		for(int j=0;j<q;++j){
			double e = 0;
			for(int k=0;k<m;++k){
				e = e + A[i][k]*B.get_entry(k,j);
			}
			//cout<<e<<endl;
			C.set_entry(i,j,e);
		}
	}
	//cout<<"leaving rh_multiply"<<endl;
	return C;
}

Matrix Matrix::operator-()const{
	Matrix B(n,m);
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			B.set_entry(i,j,-A[i][j]);
		}
	}
	return B;
}

Matrix Matrix::operator+(const Matrix& B)const{
	if(n==B.get_height() && m==B.get_width()){
		Matrix C(n,m);
		for(int i=0;i<n;++i){
			for(int j=0;j<m;++j){
				double e = A[i][j]+B.get_entry(i,j);
				C.set_entry(i,j,e);
			}
		}
		return C;
	}else{
		cerr<<"Addition Error: Dimensions do not match!"<<endl;
		Matrix C(1,1);
		return C;
	}
}

Matrix Matrix::operator-(const Matrix& B)const{
	//cout<<"in A-B"<<endl;
	Matrix C = *this+(-B);
	//cout<<"leaving A-B"<<endl;
 	return C;
}

//TO DO: ADD INDUCED MATRIX NORM
//returns 2 norm of row or column vector
double Matrix::norm_2()const{
	double ans;
	if(n==1){
		double sum = 0;
		for(int j=0;j<m;++j){
			sum = sum + A[0][j]*A[0][j];
		}
		ans = sqrt(sum);
	}else if(m==1){
		double sum = 0;
		for(int i=0;i<n;++i){
			sum = sum + A[i][0]*A[i][0];
		}
		ans = sqrt(sum);
	}else{
		//non vector objects will fail
		ans = -1;
	}
	return ans;
}

//takes the dot product of 2 column vectors
double Matrix::dot(const Matrix& v)const{
	//TO DO: ADD DIMENSION CHECKS and ADD ROW VECTORS
	double ans = 0;
	for(int i=0;i<n;++i){
		ans = ans + A[i][0]*v.get_entry(i,0);
	}
	return ans;
}

Matrix Matrix::scalar_multiply(double a)const{
	Matrix v(n,m);
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			double e = A[i][j]*a;
			v.set_entry(i,j,e);
		}
	}
	return v;
}

/*/////////////////////////////////////////////////////////////////////////
			Jacobi Operations
/////////////////////////////////////////////////////////////////////////*/

Matrix Matrix::get_Dinv()const{
	Matrix Dinv = Matrix(n,m);
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			if(i==j){
				double e = 1/A[i][j];
				Dinv.set_entry(i,j,e);
			}else{
				Dinv.set_entry(i,j,0);
			}
		}
	}
	return Dinv;
}

Matrix Matrix::get_LplusU()const{
	Matrix LplusU = Matrix(n,m);
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			if(i==j){
				LplusU.set_entry(i,j,0);
			}else{
				double e = A[i][j];
				LplusU.set_entry(i,j,e);
			}
		}
	}
	return LplusU;
}

/*/////////////////////////////////////////////////////////////////////////
			GMRES Operations
/////////////////////////////////////////////////////////////////////////*/

Matrix Matrix::get_subMatrix(int p, int q)const{
	Matrix C(p,q);	
	if(p<=n && q<=m){
		for(int i=0;i<p;++i){
			for(int j=0;j<q;++j){
				C.set_entry(i,j,A[i][j]);
			}
		}
	}else{
		cerr<<"Dimensions larger than original matrix!"<<endl;
	}
	return C;
}






