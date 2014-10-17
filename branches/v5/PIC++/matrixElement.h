#ifndef _MATRIX_ELEMENT_H_
#define _MATRIX_ELEMENT_H_

struct MatrixElement
{
	double value;
	int i;
	int j;
	int k;

	MatrixElement();
	MatrixElement(double v, int iv, int jv, int kv);

	bool equalsIndex(MatrixElement& element);
};

#endif