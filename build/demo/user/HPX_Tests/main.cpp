#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cmath> //abs
#include <string>

#include "Matrix.h"


using namespace std;

double norm_2(const vector<double>& v){
	double sum = 0;
	for(int i=0;i<v.size();++i){
		sum = sum + v[i]*v[i];
	}
	double ans = sqrt(sum);
	return ans;
}

double dot(const vector<double>& u, const vector<double>& v){
	double x = 0;
	if(u.size() == v.size()){
		for(int i=0;i<u.size();++i){
			x += u[i]*v[i];
		}
	}else{
		cerr<<"Vector lengths do not match!"<<endl;
		x = -1;
	}
	return x;
}

vector<double> v_diff(const vector<double>& u, const vector<double>& v){
	vector<double> w;	
	if(u.size() == v.size()){
		for(int i=0;i<v.size();++i){
			w.push_back(u[i]-v[i]);
		}
	}else{
		cerr<<"Vector lengths do not match!"<<endl;	
	}
	return w;
}

vector<double> sc_mult(double a, const vector<double>& v){
	vector<double> w;	
	for(int i=0;i<v.size();++i){
		w.push_back(a*v[i]);
	}
	return w;
}

void v_print(vector<double> v){
	int n = v.size();
	cout<<"[";
	for(int i=0;i<n-1;++i){
		cout<<v[i]<<", ";
	}
	cout<<v[n-1]<<"]"<<endl;
}

vector<double> foreward_eliminate(const Matrix& L, const vector<double>& b){
	int n = b.size();		
	vector<double> y(n,0);
	double t;
	for(int i=0;i<n;++i){
		t = b[i];
		for(int j=0;j<i;++j){
			t -= L.get_entry(i,j)*y[j];
		}
		y[i] = t / L.get_entry(i,i);
	}
	return y;
}

vector<double> back_substitute(const Matrix& U, const vector<double>& y){
	int n = y.size();
	vector<double> x(n,0);
	double t;
	for(int i=n-1;i>=0;--i){
		t = y[i];
		for(int j=i+1;j<n;++j){
			t -= U.get_entry(i,j)*x[j];
		}
		x[i] = t / U.get_entry(i,i);
	}
	return x;
}

void gauss(const Matrix& A, const vector<double>& b, vector<double>& xa){
	int n = A.get_height();	
	//factor
	Matrix L(n,n);
	Matrix U(A);
	double e;
	for(int j=0;j<n-1;++j){
		for(int i=j+1;i<n;++i){//
			e = U.get_entry(i,j) / U.get_entry(j,j);
			L.set_entry(i,j,e);
			for(int k=0;k<n;++k){
				e = U.get_entry(i,k)-L.get_entry(i,j)*U.get_entry(j,k);
				U.set_entry(i,k,e);
			}
		}
		L.set_entry(j,j,1);
	}
	L.set_entry(n-1,n-1,1);
	
	//Solve LUx=b
	xa = back_substitute(U,foreward_eliminate(L,b)); 
}

//Reference for Pseudocode (MATLAB / GNU Octave):
//https://en.wikipedia.org/wiki/Generalized_minimal_residual_method

//Input empty vectors q,h
void arnoldi(Matrix A, Matrix Q, int k, vector<double>& q, vector<double>& h){
	//Q.print();
	q = A.mv_mult(Q.get_col(k));
	for(int j=0;j<k+1;++j){ //changed k->k+1
		h.push_back(dot(q,Q.get_col(j)));//h[k]
		for(int i=0;i<q.size();++i){
			q[i] = q[i] - h[j]*Q.get_entry(i,j);
		}
	}
	h.push_back(norm_2(q));//h[k+1]
	//normalize q
	q = sc_mult(1/(h[k+1]),q);
}

void givens_rotation(double v1, double v2, double& cs, double& sn){
	if(v1 == 0){
		cs = 0;
		sn = 1;
	}else{
		double t = sqrt(v1*v1+v2*v2);
		cs = abs(v1) / t;
		sn = cs * v2 / v1;
	}
}

void apply_givens_rotation(vector<double>& h, vector<double>& cs, vector<double>& sn, int k){
	//apply for ith col
	for(int i=0;i<k-1;++i){
		double temp = cs[i]*h[i] + sn[i]*h[i+1];
		h[i+1] = -sn[i]*h[i] + cs[i]*h[i+1];
		h[i] = temp;
	}
	
	//update the next sin cos values for rotation
	givens_rotation(h[k], h[k+1], cs[k], sn[k]);

	//eliminate H(k+1,k)
	h[k] = cs[k]*h[k] + sn[k]*h[k+1];
	h[k+1] = 0;
}



//Only works for A:nxn
void gmres(const Matrix& A, const vector<double>& b, vector<double>& xa,  double tol, int maxIter){
	int n = A.get_height();
	int m = maxIter;
	
	//use xa as the initial vector
	vector<double> r = v_diff(b,A.mv_mult(xa));
	double b_norm = norm_2(b);
	double r_norm = norm_2(r);
  	double error = r_norm/b_norm;

	//initialize the 1D vectors
	vector<double> sn(m,0);
	vector<double> cs(m,0);
	vector<double> e1(n,0);
	e1[0] = 1;
	//e = [error] ?
	//Allocate matrices to maximum size
	Matrix Q(n,m+1);
	Matrix H(m+1,m);
	
	
	Q.set_col(0,sc_mult(1/r_norm,r));
	vector<double> Beta = sc_mult(r_norm,e1);

	int totalIter;
	for(int k=0;k<m;++k){
		//run arnoldi
		vector<double> q,h;
		arnoldi(A, Q, k, q, h); //ERROR HERE
		H.set_col(k,h);
		Q.set_col(k+1,q);
		
		//H.print(); Q.print();

		//eliminate the last element in H ith row and update rotation
		apply_givens_rotation(h, cs, sn, k);
		H.set_col(k,h);
		
		//update the residual
		Beta[k+1] = -sn[k]*Beta[k];
		Beta[k] = cs[k]*Beta[k];
		error = abs(Beta[k+1]) / b_norm;
	
		//save the error?
	
		if( error <= tol){
			totalIter = k+1; //INDEX ERROR?
			cout<<"Took "<<totalIter<<" Iterations"<<endl;
			break;
		} 
	}
	
	//calculate the result:

	//truncate the system
	Matrix H_sub = H.get_subMatrix(totalIter,totalIter);
	Beta.resize(totalIter);

	//H_sub.print(); v_print(Beta);	

	//solve upper triangular
	vector<double> y = back_substitute(H_sub, Beta);
	//update the solution
	Matrix Q_sub = Q.get_subMatrix(n,totalIter);

	//Q_sub.print(); v_print(y);

	vector<double> Qy = Q_sub.mv_mult(y);
	for(int i=0;i<n;++i){
		xa[i] += Qy[i];
	}
}

//end reference

void jacobi(const Matrix& A, const Matrix& b, Matrix& xa, double tol, int maxIter){
	//cout<<"Begin Jacobi"<<endl;
	//form Jacobi iteration using B,d
	Matrix Dinv = A.get_Dinv();	 	//Dinv.print();
	Matrix negDinv = -Dinv;		 	//negDinv.print();
	Matrix LplusU = A.get_LplusU();  	//LplusU.print();

	Matrix B = negDinv.rh_multiply(LplusU); //B.print();
	Matrix d = Dinv.rh_multiply(b);
	Matrix r;
	for(int k=0;k<maxIter;++k){
		//check residual
		r = A.rh_multiply(xa)-b;	//r.print();
		if(r.norm_2()<tol){
			cout<<"Took "<<k<<" Iterations"<<endl;
			return;
		}
		
		xa = B.rh_multiply(xa) + d;

		//xa.print();

	}
	//return 0 matrix if maxIter exceed
	cerr<<"Maximum Iterations Exceeded!"<<endl;
}

Matrix conjugateGradient(const Matrix& A, const Matrix& b, Matrix& xa, double tol, int maxIter){
	//cout<<"Begin ConjGrad"<<endl;
	Matrix r = b-A.rh_multiply(xa);
	Matrix d(r); //direction
	double a; Matrix r_next(b.get_height(),1); double B;
	for(int k=0;k<maxIter;++k){
		//check residual
		if(r.norm_2()<tol){
			cout<<"Took "<<k<<" Iterations"<<endl;
			return xa;
		}
		
		a = r.dot(r) / d.dot(A.rh_multiply(d)); //step length
		xa = xa + d.scalar_multiply(a); //xa.print();
		r_next = r - A.rh_multiply(d).scalar_multiply(a); 
		B = r_next.dot(r_next) / r.dot(r);
		d = r_next + d.scalar_multiply(B);
		
		Matrix tmp(r_next);
		r = tmp;
	}
	//return 0 matrix if maxIter exceed
	cerr<<"Maximum Iterations Exceeded!"<<endl;
	return Matrix(1,1);
}

Matrix constructMatrix(const vector<double>& s,const vector<double>& t){
	//TO DO: CASE WHERE S AND T ARE DIFFERENT SIZES!
	int n = s.size(); int m = t.size();
	Matrix A(n,m);
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			//perform a pseudo-ish inverse
			double e; 
			if(abs(s[i]-t[j])<.0000001){e = n;}
			else{e = 1/abs(s[i]-t[j]);}
			A.set_entry(i,j,e);
		}
	}
	return A;
}

int main(int argc, char* argv[]){
	/*//initialize hpx
	if (hpx_init(&argc, &argv)) {
    		hpx_print_help(); 
    		return -1;
  	}*/
	
	if(argc != 4){ 
		//FIX: WRITE BETTER ERROR MESAGE
		cerr<<"Usage: string IterationType double tol int maxIter 			"<<endl;
		return -1;
	}
	//check for valid IterationType
	int itType;
	string str(argv[1]);
	if(str == (string)"Gauss"){itType = 0;}
	else if (str == (string)"Jacobi"){itType = 1;}
	else if (str == (string)"ConjugateGradient"){itType = 2;}
	else if (str == (string)"GMRES"){itType = 3;}
	else{
		cerr<<"Invalid Algorithm:"<<endl;
		cerr<<"Try: Gauss, Jacobi, ConjugateGradient, GMRES"<<endl;
		return -1;
	}

	double tol = atof(argv[2]); int maxIter = atoi(argv[3]);

	//create test A and b
	int n = 10;
	//Matrix A(n,n);
	vector<double> s;
	vector<double> t;
	for(int i=0;i<n;++i){
		cout<<i<<" ";
		s.push_back(i);
		t.push_back(i);
	}
	cout<<endl<<endl;

	Matrix A = constructMatrix(s,t);	A.print();

						

	//create soln vector
	Matrix xe(n,1,1); 			xe.print();
	Matrix b = A.rh_multiply(xe);		b.print();

	//Perform Iterative Algorithm
	Matrix xa(n,1,0);/*guess 0 vector*/	xa.print();
	//
	vector<double> b_vect = b.toVector();	
	vector<double> xa_vect = xa.toVector();

	switch(itType){
		case 0: gauss(A,b_vect,xa_vect);
			break;
		case 1: jacobi(A,b,xa,tol,maxIter);
			break;
		case 2: conjugateGradient(A,b,xa,tol,maxIter);
			break;
		case 3: gmres(A,b_vect,xa_vect,tol,maxIter);
			break;
	}
	if(itType == 1 || itType == 2){xa.print();}
	else{v_print(xa_vect);}
	
	
	
	cout<<"Hello World"<<endl;
	return 0;
}

