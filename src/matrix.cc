
#include "matrix.h"

// default constructor
matrix::matrix() {
	value = new double*[row=1];
	value[0] = new double[col=1];
	value[0][0] = 0.0;
}	

// constructor
matrix::matrix(const int idx1, const int idx2) {         
#ifdef MATRIXBOUNDS
	if(idx1<1 || idx2<1) 
		error("[Matrix] Error: matrix/vector size must exceed 0.\n");
#endif
	row = idx1;
	col = idx2;
	value = new double*[row];
	for (int i=0;i<row;i++){
		value[i]=new double[col];
		for (int j=0;j<col;j++) 
			value[i][j]=0.0;
	}
//	matrixConst += 1;
//  	Rprintf("mat: %d \n",matrixConst);
}

//matrix to matrix copy constructor
matrix::matrix(const matrix &a) {
	value = new double*[row=a.row];
	for (int i=0;i<row;i++) {
		value[i]=new double[col=a.col];
		for (int j=0;j<col;j++)
			value[i][j]=a.value[i][j];
	}
	/* 
	 * matrixConst += 1;
	 * Rprintf("mat: %d \n",matrixConst);
	 */
}

//destructor
matrix::~matrix() {  
	for (int i=0;i<row;i++) 
		delete [] value[i];
	delete [] value;
}

void matrix::reset(const int idx1, const int idx2) {
#ifdef MATRIXBOUNDS
	if (idx1<1 || idx2<1) 
		error("[Matrix] Error: reset matrix/vector size must exceed 0.\n");
#endif
	for (int i=0;i<row;i++) 
		delete [] value[i];
	delete [] value;
	value = new double*[row=idx1];
	for (int i=0;i<row;i++) {
		value[i]=new double[col=idx2];
		for (int j=0;j<col;j++) 
			value[i][j]=0.0;
	}
}

//matrix access
double& matrix::operator() (const int x,const int y) {
#ifdef MATRIXBOUNDS	
	if(x<1 || x>row || y<1 || y> col) error("[Matrix] Error: matrix out of range.\n");
#endif
	return(value[x-1][y-1]);
}
		  
//vector access
double& matrix::operator() (const int x) {
#ifdef MATRIXBOUNDS			
	if(row==1) {
		if(x<1 || x>col) error("[Matrix] Error: rowvector out of range.\n");
	}
	else {
		 if(col==1) {
		 	if (x<1 || x>row) error("[Matrix] Error: colvector out of range.\n");
		 }
		 else error("[Matrix] Error: matrix adressed as vector.\n");
	}
#endif
	if (row==1) return(value[0][x-1]);
	else {
		if (col==1) return(value[x-1][0]);
	}
}

//cast
matrix::operator double() {
#ifdef MATRIXBOUNDS	
	if(!(col==1 && row==1)) error("[Matrix] Error: incorrect matrix to double cast.\n");
#endif
	return value[0][0];
}

// = (assignment)
matrix& matrix::operator=(const matrix &a) {
	double **temp;
	temp = new double*[a.row];
	for (int i=0;i<a.row;i++) {
		temp[i]=new double[a.col];
		for (int j=0;j<a.col;j++) 
			temp[i][j]=a.value[i][j];
	}
	for (int i=0;i<row;i++) 
		delete [] value[i];
	delete [] value;
	row=a.row;
	col=a.col;
	value=temp;
	return *this;
}

// = assignment from double
matrix& matrix::operator=(const double a) {
#ifdef MATRIXBOUNDS
	if (row!=1 || col!=1) error("[Matrix] Error: incorrect scalar assigment to matrix/vector.\n");
#endif
	value[0][0]=a;
	return *this;
}

// row and col operations

//returns a colvector of rowsums
matrix matrix::rowsums() {			
	matrix target(row);
	for(int i=0;i<row;i++) {
		for(int j=0; j<col; j++) {
			target.value[i][0] += value[i][j];
		}
	}
	return(target);
}

//returns a rowvector of colsums
matrix matrix::colsums(){			
	matrix target(col);
	for(int j=0;j<col;j++) {
		for(int i=0; i< row; i++) {
			target.value[0][j] += value[i][j];
		}
	}
	return(target);
}

double matrix::msum() { 
	if(row==1 && col==1) 
		return(value[0][0]);	
	double csum=0.0;
	if(col == 1) {
		for(int i=0;i<row;i++)
			csum += value[i][0];
		return(csum);
	}
	if(row==1)  {
		for(int i=0;i<col;i++)
			csum += value[0][i];
		return(csum);
	}
	else 
		error("[Matrix] sum only defined for row or col vector.\n");
}

matrix matrix::rown(int rowNr) {
	matrix target(1,col);
	for(int j=0;j<col;j++)
		target.value[0][j] = value[rowNr-1][j];
	return(target);
}

matrix matrix::coln(int colNr) {
	matrix target(row,1);
	for(int j=0;j<row;j++)
		target.value[j][0] = value[j][colNr-1];
	return(target);
}

// +
matrix matrix::operator+(const matrix &b) {
#ifdef MATRIXBOUNDS
	if (!(row==b.row && col==b.col)) error("[Matrix] Error: incompatible matrices + .\n");
#endif
	matrix target(row,col);
	for (int i=0;i<row;i++) 
		for (int j=0;j<col;j++) 
		target.value[i][j]=value[i][j]+b.value[i][j];
	return(target);
}

matrix matrix::operator+(const double b) { 
#ifdef MATRIXBOUNDS
    if (col!=1 || row!=1) error("[Matrix] Error: cannot add scalar to matrix/vector.\n");
#endif
	matrix target(1,1);
	target.value[0][0]=value[0][0]+b;
	return(target);
}

// -
matrix matrix::operator-(const matrix &b) {
#ifdef MATRIXBOUNDS
	if (!(row==b.row && col==b.col)) error("[Matrix] Error: incompatible matrices - .\n");
#endif
	matrix target(row,col);
	for (int i=0;i<row;i++) 
		for (int j=0;j<col;j++) 
			target.value[i][j]=value[i][j]-b.value[i][j];
	return(target);
}

matrix matrix::operator-(const double b) {
#ifdef MATRIXBOUNDS
	if (col!=1 || row!=1) error("[Matrix] Error: cannot substract scalar from matrix/vector.\n");
#endif
	matrix target(1,1);
	target.value[0][0]=value[0][0]-b;
	return(target);
}

// *  matrix products
matrix matrix::operator*(const matrix &b) {
	if (row==1 && col==1) { //left hand scalar
		matrix target(b.row,b.col);
		for (int i=0;i<b.row;i++) for (int j=0;j<b.col;j++) 
			target.value[i][j]=b.value[i][j]*value[0][0];
		return(target);
	}
	else { //right hand scalar
		if (b.row==1 && b.col==1) {
			matrix target(row,col);
			for (int i=0;i<row;i++) 
				for (int j=0;j<col;j++) 
					target.value[i][j]=value[i][j]*b.value[0][0];
			return(target);
		}
		else { //ordinary matrix muliplication
			if (col==b.row) {
				matrix target(row,b.col);
				for (int i=0;i<row;i++) 
					for (int j=0;j<b.col;j++) 
						for (int k=0;k<col;k++)
						target.value[i][j]=target.value[i][j]+value[i][k]*b.value[k][j];
					return(target);
			}
			else error("[Matrix] Error: incompatible matrices *.\n");
		}
	}
}

matrix matrix::operator*(const double b) {
	matrix target(row,col);
	for (int i=0;i<row;i++) for (int j=0;j<col;j++) 
		target.value[i][j]=value[i][j]*b;
	return(target);
}

// had hadamard product

matrix had(const matrix &a,const matrix &b) {
#ifdef MATRIXBOUNDS
	if(a.row!=b.row || a.col!=b.col) error("[Matrix] Error: nonconformable matrices in hadamard matrix product.\n");
#endif
	matrix target(a.row,a.col);
	for (int i=0;i<a.row;i++) 
		for (int j=0;j<a.col;j++) 
			target.value[i][j]=a.value[i][j]*b.value[i][j];
	return(target);
}

// / (division)			 
matrix matrix::operator/(const matrix &b) { // returns a/b b is scalar
#ifdef MATRIXBOUNDS
	if (b.col!=1 || b.row!=1) error("[Matrix] Error: matrix incorrectly adressed as scalar in division.\n");
	if (b.value[0][0]==0.0) error("[Matrix] Error: division by zero.\n");
#endif
	matrix target(row,col);
	for (int i=0;i<row;i++) 
		for (int j=0;j<col;j++) 
			target.value[i][j]=value[i][j]/b.value[0][0];
	return(target);
}

matrix matrix::operator/(const double b) { // returns a/b 
#ifdef MATRIXBOUNDS
	if (b==0.0) error("[Matrix] Error: division by zero.\n");
#endif
	matrix target(row,col);
	for (int i=0;i<row;i++) 
		for (int j=0;j<col;j++) 
			target.value[i][j]=value[i][j]/b;
	return(target);
}

//transpose
matrix transpose(const matrix &a) {
	matrix target(a.col,a.row);
	for (int i=0;i<a.row;i++) 
		for (int j=0;j<a.col;j++) 
			target.value[j][i]=a.value[i][j];
	return(target);
}

double max(matrix a) {
	if(!(a.row==1||a.col==1)) error("[Matrix] max only defined for row or col vector.\n");
	int maxidx=1;
	int idx=0; 
	double max=a(1);
	if(a.row==1) idx=a.col;
	else idx=a.row;
	for(int i=1; i<=idx; i++) {
		if(a(i)>max) max = a(i);
	}
	return(max);
}

int argmax(matrix a) {
	if(!(a.row==1||a.col==1)) error("[Matrix] argmax only defined for row or col vector.\n");
	double max=a(1);
	int maxidx=1;
	int idx=0; 
	if(a.row==1) idx=a.col;
	else idx=a.row;
	for(int i=1; i<=idx; i++) {
		if(a(i)>max) {
			max=a(i);
			maxidx=i;
		}
	}
	return(maxidx);
}

void matrix::normalize(void) {
	if(!(row==1||col==1)) error("[Matrix] normalize only defined for row or col vector.\n");
	double ms=msum();
	if(row>1) {
		for(int i=0; i<row; i++) value[i][0] /= ms;
	}
	else for(int i=0; i<col; i++) value[i][0] /= ms;
}

void matrix::print(void) {
	for(int i=0;i<row; i++) {
		for(int j=0;j<col; j++) {
			Rprintf(" %f",value[i][j]); 
		}
		Rprintf("\n");
	}
	Rprintf("\n");
}
