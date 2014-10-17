#include "matrixElement.h"

MatrixElement::MatrixElement() {
	value = 0;
	i = 0;
	j = 0;
	k = 0;
}

MatrixElement::MatrixElement(double v, int iv, int jv, int kv) {
	value = v;
	i = iv;
	j = jv;
	k = kv;
}

bool MatrixElement::equalsIndex(MatrixElement& element) {
	if(i != element.i) return false;
	if(j != element.j) return false;
	if(k != element.k) return false;
	
	return true;
}