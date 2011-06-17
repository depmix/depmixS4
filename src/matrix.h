
#ifndef MATRIX
#define MATRIX 1

#include <stdio.h>
#include <stdlib.h>

extern "C" {
	#include <R.h>
	static int matrixConst; 
}

#define MATRIXBOUNDS 1

class matrix {
	public:
		int row,col;
		double  **value;
 		
		
	public:
		// constructors - deconstructors
		matrix();          								// default constructor
		matrix(const int idx1,const int idx2=1);		// constructor
		matrix(const matrix &a);						// copy constructor
		~matrix();										// destructor
		
		void reset(const int idx1, const int idx2=1);
		
		// accesss & cast
		double &operator()(const int x,const int y);	// matrix access
		double &operator()(const int x); 				// vector access
		operator double();								// cast matrix to double
		
		// =
		matrix& operator=(const matrix &a);				// assignment operator
		matrix& operator=(const double a);				// assignment operator
		
		// row & col operations
		friend inline int rows(const matrix &a) {return a.row;} // returns # of rows
		friend inline int cols(const matrix &a) {return a.col;} // returns # of colums

		matrix rowsums();			//returns a colvector of rowsums
		matrix colsums();			//returns a rowvector of colsums
		double msum();				//returns the sum of a col or row vector
		
 		matrix rown(int rowNr);		// returns a rowvector rowNr 
		matrix coln(int colNr);		// returns a colvector colNr 
		
		// +
		matrix operator+(const matrix &b); 	// returns a+b
		matrix operator+(const double b);  	// if a is scalar
		
		// -
		matrix operator-(const matrix &b);	// returns a-b
		matrix operator-(const double b); 	// if a is scalar
		
		// *	matrix product
		matrix operator*(const matrix &b); 	// returns a*b
		matrix operator*(const double b);  	// b scalar
		
		// had hadamard matrix products
		friend matrix had(const matrix &a,const matrix &b);	// returns hadamard product of a and b
		
		// / (division)
		matrix operator/(const matrix &b); 			// returns a/b b is scalar
		matrix operator/(const double b);  			// returns a/b
		
		// other operations/stuff
		friend matrix transpose(const matrix &a); 	// returns transpose
		friend double max(matrix a);				//returns the largest element of a col or row vector
		friend int argmax(matrix a);				//returns the index of the above
		void normalize(void);						//normalize a row or col vector such that its values represent probs
		
		// i/o 
		void print(void);							// screen output
		
};

#endif

