#include <iostream>
#include <sstream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

const mat NaturalSpline {{-1.0f/6.0f,	   0.5f,  	 -0.5f, 1.0f/6.0f},
						 {	    0.5f, 	  -1.0f,  	  0.5f,	     0.0f},
					  	 {	   -0.5f,  	   0.0f,      0.5f,		 0.0f},
						 { 1.0f/6.0f, 4.0f/6.0f, 1.0f/6.0f,		 0.0f}};
								 
const mat CubicMiddles {{-1.0f,  3.0f, -3.0f,  1.0f},
                        { 2.5f, -5.5f,  3.5f, -0.5f},
                        {-2.0f,  2.0f,  0.0f,  0.0f},
                        { 0.5f,  0.5f,  0.0f,  0.0f}};
                              
const mat Middles_3 {{ 0.5f, -1.0f, 0.5f},
				     {-1.0f,  1.0f, 0.0f},
					 { 0.5f,  0.5f, 0.0f}};

int Pascal(int i, int n)
{
	if ((i<0)||(i>n)) return 0;
	else if ((i==0)||(i==n)) return 1;
	else return Pascal(i-1, n-1) + Pascal(i, n-1);
}

mat BezierMatrix(int n) {
	vector<int> Values;
	for(int j=0; j<n; j++) {
		for(int i=0; i<=j; i++) {
			Values.push_back(pow(-1,n-1-j)*Pascal(j,n-1)*Pascal(i,j));
		}
	}
	mat M = zeros(n,n);
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			if ((i+j)>=n) M(i,j)=0;
			else M(i,j)=Values[(((i+j)*(i+j+1))/2)+j];
		}
	}
	return M;
}

int main()
{
	mat toFile;
	ostringstream fileName;
	
	// Les matrices de Hermite
	for (int j=0; j<=24; j++) {
		float alpha = (float)j/12.0f;
		toFile = {{		2.0f-alpha,		  alpha, -alpha-2.0f,  alpha},
				  {2.0f*alpha-3.0f, -2.0f*alpha,  alpha+3.0f, -alpha},
				  {			-alpha,		  alpha, 		0.0f,   0.0f},
				  {			  1.0f,		   0.0f, 		0.0f,   0.0f}};
		fileName << "./Matrices/Hermite/Hermite_" << j << ".mat";
		toFile.save(fileName.str());
		fileName.str("");
		fileName.clear();
		toFile.reset();
	}
	
	// Les matrices de Bezier 
	for(int k=1; k<=25; k++) {
		toFile = BezierMatrix(k);
		fileName << "./Matrices/Bezier/Bezier_" << k << ".mat";
		toFile.save(fileName.str());
		fileName.str("");
		fileName.clear();
		toFile.reset();
	}
	
	// Les matrices de CatmullRom
	for(int i=0; i<=24; i++) {
		float alpha = (float)i/12.0f;

		toFile = {{    -alpha, 2.0f-alpha,      alpha-2.0f,  alpha},
				  {2.0f*alpha, alpha-3.0f, 3.0f-2.0f*alpha, -alpha},
				  {    -alpha,       0.0f,           alpha,   0.0f},
				  {      0.0f,       1.0f,            0.0f,   0.0f}};
				  
		fileName << "./Matrices/CatmullRom/CatmullRom_" << i << ".mat";
		toFile.save(fileName.str());
		fileName.str("");
		fileName.clear();
		toFile.reset();
	}
	
	// La matrice des BSplines cubiques non-rationnelles et uniformes
	fileName << "./Matrices/NaturalSpline.mat";
	NaturalSpline.save(fileName.str());
	fileName.str("");
	fileName.clear();
		
	// La matrice des milieux quadratiques
	fileName << "./Matrices/Middles_3.mat";
	Middles_3.save(fileName.str());
	fileName.str("");
	fileName.clear();
		
	// La matrice des milieux cubiques 3 milieux
	fileName << "./Matrices/Middles_4.mat";
	CubicMiddles.save(fileName.str());
	fileName.str("");
	fileName.clear();
	
	return 0;
}