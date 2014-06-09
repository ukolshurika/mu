#include <iostream>
#include <fstream>

#include "membrane.h"

using namespace std;

int main(){
  double h0 = 2/100.0;
  double q = 2.65*1000/88.3/1000000;
  double n = 3.4;
  double mu = 0.1;

  Membrane m(q, h0, n, mu);
  ofstream free_data("data/free_new.dat");

  m.free(999);
  m.constrained(1000);


  for(auto i = m.t_free_.cbegin(); i != m.t_free_.cend(); ++i){
    free_data<< i -> first << ' '<< i->second << endl;
  }


  cerr << "C0MPLETE" << endl;
  return 0;
}
