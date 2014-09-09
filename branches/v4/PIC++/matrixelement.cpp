#include "matrixElement.h"

MatrixElement::MatrixElement(double v, int iv, int jv, int kv, int lv) {
	value = v;
	i = iv;
	j = jv;
	k = kv;
	l = lv;
}

bool MatrixElement::equalsIndex(MatrixElement& element) {
	if(i != element.i) return false;
	if(j != element.j) return false;
	if(k != element.k) return false;
	
	return true;
}