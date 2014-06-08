#include "membrane.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

#include "bound.h"
#include "matrix.h"
#include "simpson.h"
#include "utils.h"

using namespace std;

namespace{
  const double kSqrt3 = sqrt(3);
  const int kSimpsonStep = 99;
}

struct SFunctor{
  SFunctor(){};
  double operator()(double x) const{
    return sqrt(1+Matrix::kKBSquare*x*x);
  }
};

struct Free {
  Free(double h0, double q, double n):q_(q), h0_(h0), n_(n){};
  double operator()(double alpha) const{
    return (1/alpha-1/tan(alpha))*pow(((2*h0_*sin(alpha)*sin(alpha))/(kSqrt3*q_*alpha)-1), n_);
  }

  private:
  double q_;
  double h0_;
  double n_;
};

Membrane::Membrane(double q, double h0, double n, double mu):q_(q), h0_(h0), n_(n), mu_(mu){
  alpha1_ = 0.41; // from Maxima flexible step(boolsh it just get it from terraud)
  alpha2_ = Bound((*this)).Alpha(Matrix::RZero());
  cerr << alpha2_ << endl;
  h1_ = sin(alpha2_)/alpha2_*h0_;
}

void Membrane::free(int steps){
  double dalpha = (alpha2_ - alpha1_)/steps;
  double t;
  Free f(h0_, q_, n_);

  vector<pair<double, double>> v;

  for(double a = alpha1_; a < alpha2_; a+=dalpha){
    t = Simpson::Integrate(a, a+dalpha, kSimpsonStep, f);
    v.push_back(make_pair(t, a));
  }

  t_free_.clear();
  double offset = 0;

  for (auto it = v.begin(); it != v.end(); ++it) {
    t_free_.push_back(make_pair((it->first + offset), it->second));
    offset += it->first;
  }
}

void Membrane::prepare_x_k(){
  cerr << "START PREPARE\n";
  for(double x=1; x>0; x-=0.0001){
    s_x_.push_back(make_pair(Simpson::Integrate(x, 1, kSimpsonStep, SFunctor()), x));
    // cout << s_x_.back().first << " " << x << endl;
  }
  sort(s_x_.begin(),s_x_.end());
  cerr << "END PREPARE\n";
}

double Membrane::find_x(double s){
  //cerr << s << endl;
  auto it = s_x_.begin();
  for(it; it-> first < s; ++it){};
  double s_1 = it->first;
  double x_1 = it->second;
  //cerr << s_1 << '!' << x_1 << endl;
  ++it;
  double s_2 = it->first;
  double x_2 = it->second;
  //cerr << s_2 << '!' << x_2 << "!!" << (s-s_1)/(s_2-s_1)*(x_2-x_1)+x_1<< endl;
  return (s-s_1)/(s_2-s_1)*(x_2-x_1)+x_1;
}

void Membrane::constrained(int steps){
  prepare_x_k();

  ofstream h_data("data/h.dat");
  ofstream sigma_data ("data/s.dat");
  ofstream x_data ("data/x.dat");


  double init_sigma_k = q_/h1_; // Q*RHO/H (from free stadia)
  double dt = 1000000;
  double x, h=0;

  vector <double> sigma_k(steps, init_sigma_k), sigma_k1(steps, 0.0),
                  ds_k(steps, 0.0), ds_k1(steps, 0.0), 
                  delta_ds_k(steps, 0.0), delta_ds_k1(steps, 0.0),
                  h_k(steps, h1_/h0_), h_k1(steps, 0.0),
                  d_rho_k(steps, 0.0), rho_k(steps, 0.0), 
                  alpha_k(steps, alpha2_), d_alpha_k(steps, 0),
                  sk(steps, 0), sk1(steps, 0);
  for(auto t = 1; t < steps ; ++t){
    
    for(auto i = 1; i < t; ++i){
      delta_ds_k1[i] = pow((sigma_k[i-1]+sigma_k[i])/(4/sqrt(3)-(sigma_k[i-1]+sigma_k[i] )), n_)*ds_k[i]*dt;
      // cerr << sigma_k[i-1] << endl;
      DCHECK(delta_ds_k1[i]>=0);
    }

    for(auto i = 0; i< t; ++i){
      ds_k1[i] = ds_k[i] + delta_ds_k1[i];
      DCHECK(ds_k1[i]>=0);
    }
    
    fill(sk.begin(), sk.end(), 0);
    for(auto i = 0; i<t; ++i){ 
      for(auto j = 0; j<=i; ++j) 
        sk[i] += ds_k1[j];
        DCHECK(sk[i]>=0);
    }

    double x_prev = 1.0, dx;
    for(auto i = 0; i<t; ++i){
      x = Membrane::find_x(sk[i]);
      dx = x-x_prev;
      if(dx < 1e-8) dx=0;
      DCHECK(x <= x_prev);
      rho_k[i] = sqrt(x*x+x*x/Matrix::kKBSquare/pow(x, 2*Matrix::kK-2));
      d_rho_k[i] = (x*x+(2-Matrix::kK)*pow(x, 3-2*Matrix::kK)/Matrix::kKBSquare)/rho_k[i]*dx;
      alpha_k[i] = M_PI_2 - atan(1/(Matrix::kK*Matrix::kB*pow(x, Matrix::kK-1)));
      double k = Matrix::kK;
      double b = Matrix::kB;
      d_alpha_k[i] =   (k*(k-1)*b*pow(x, k-2)/(1+Matrix::kKBSquare*pow(x, 2*k-2)))*dx;
      x_prev = x;
    }
    
    ofstream x_data ("data/x"+to_string(t)+".dat");

    for(auto i = 0; i<t; ++i){
      h_k1[i] = h_k[i]-h_k[i]*pow(1/(1-sqrt(3)/2*q_*rho_k[t-1]/h0_/h_k[t-1])-1, n_)*dt;
      cerr  << h_k[i]*pow(1/(1-sqrt(3)/2*q_*rho_k[t-1]/h0_/h_k[t-1])-1, n_)*dt<< endl;
      DCHECK(h_k1[i]>0);
      x_data << i << " "<< h_k1[i]<< endl;
    }

    ds_k1[t] = alpha_k[t-1]*rho_k[t-1]*pow(1/(1-sqrt(3)/2*q_/h0_/h_k[t-1])-1, n_)*dt -alpha_k[t-1]*d_rho_k[t-1]-rho_k[t-1]*d_alpha_k[t-1]-d_rho_k[t-1]*d_alpha_k[t-1];
    DCHECK(ds_k1[t]>0);

    x = Membrane::find_x(sk[t-1] + ds_k1[t]);    
    rho_k[t] = sqrt(x*x+x*x/Matrix::kKBSquare/pow(x, 2*Matrix::kK-2));
    alpha_k[t] = M_PI_2 - atan(1/(Matrix::kK*Matrix::kB*pow(x, Matrix::kK-1)));

    /** WARNING!!! THINK ABOUT HK+1K+1 !!!*/
    h=0;

    for(auto i = 0; i<t; ++i){
      h+=h_k1[i]*ds_k1[i];
     // cerr << h << " @"<<endl;
    }

    // // cerr << ds_k1[t] << "===="  << t<< endl;
    // if(t>1) 

      // h_k1[t]=(1-h/h0_)/(alpha_k[t]*rho_k[t]);
    // else 
      h_k1[t] = h_k1[t-1];
    // cerr << h_k1[t] << "@!" << endl;
    // DCHECK(abs(1 - h/h0_ - h_k1[t]*alpha_k[t]*rho_k[t])< 0.1);

    // cerr << h << endl;
    
    DCHECK(h_k1[t]>0);
	 
    //cerr << h_k1[t]<< "#" << t << " "<<(rho_k[t-1])<< endl;

    sigma_k1[t] = sqrt(3)/2*q_*rho_k[t]/h0_/h_k1[t-1];
    DCHECK(sigma_k1[t]>0);
    
    for(auto i = t-1; i>0; --i){
      sigma_k1[i] = sigma_k1[i+1]*h_k1[i+1]/h_k1[i]-mu_*ds_k1[i+1]*q_/h0_/h_k1[i];
      // cerr << h_k1[0]<< " " << h_k1[i] << " "<< mu_*ds_k1[i+1]*q_/h0_/h_k1[i] << endl;  
      DCHECK(sigma_k1[i]>0);
    
    }
    
    // cout << "!" << t << " "<< h_k1[t-1] << endl;
    
    h_data  << t_free_.back().first + (t)*dt << ' ' << h_k1[0] << endl;
    sigma_data << t_free_.back().first + t*dt << ' ' << sigma_k1[t]*1000/*sqrt(3)/2*q_/exp(-2/M_PI*ds_k1[t])*/<<endl;

    /* --- SWAPPING --- */
    sigma_k.swap(sigma_k1);
    ds_k.swap(ds_k1);
    delta_ds_k.swap(delta_ds_k1);
    h_k.swap(h_k1);
    // cerr << t << endl;

}
}
