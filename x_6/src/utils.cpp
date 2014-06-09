#include "utils.h"

bool util::eql(double a, double b){
	return abs(a-b) < kEpsilon;
};

bool util::IsNan(double a){
	return a != a;
};
