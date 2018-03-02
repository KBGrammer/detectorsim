#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

double cs(double x);
double is(double x);
double pe(double x);
double npp(double x);
double epp(double x);
double mu_i(double x, double array[], int nterms);
double mu_total(double x);
double cs_arr[4] = {-2.84582683378747,-1.98187129924843,-0.0528598321128271,0.0424806763283969};
double is_arr[4] = {-1.27668885411202,-0.491939420500784,-0.150907754455624,0.031031257379362};
double pe_arr[6] = {-2.35519828168804,-2.06429208226674,0.698305851759179,0.00502304809343589,
		    -0.154102421097777,-0.00918861370146097};
double npp_arr[6] = {-4.98224267460286,17.0346577242561,-45.7527308923475,68.7855291590646,
		     -52.4572355587654,15.8136280377002};
double epp_arr[11] = {-5576.265643212,90425.2750860535,-652856.065045285,2758516.73545972,
		      -7551838.66072215,13997251.3011271,-17792280.3704884,15319944.7241945,
		      -8554643.11593555,2798482.21068001,-407429.30716225};

int main()
{
  ofstream fout;
  fout.open("csfunc.txt");
  ifstream fin;
  fin.open("attenuation.txt");
  double var = 0;
  while(!fin.eof()) {
    double en,temp,temp0;
    fin >> en >> temp >> temp >> temp >> temp >> temp >> temp0 >> temp;
    double mcs = mu_i(en,cs_arr,4);
    double mis = mu_i(en,is_arr,4);
    double mpe = mu_i(en,pe_arr,6);
    double mn = mu_i(en,npp_arr,6);
    double me = mu_i(en,epp_arr,11);
    double total = mu_total(en);
    fout << en << '\t' << mcs << '\t' << mis << '\t' << mpe;
    fout << '\t' << mn << '\t' << me << '\t' << total << '\t' << (temp0-total) << endl;
    var += (total-temp0)*(total-temp0);
  }
  cout << "Sum of squares: " << var << endl;
  cout << "mu(488kev): " << mu_total(.488)/0.1004 << endl;
  cout << "mu(662kev): " << mu_total(.662)/.07754 << endl;
  cout << "mu(1000kev): " << mu_total(1.0)/.05848 << endl;
  cout << "mu(2200kev): " << mu_total(2.2)/.03997 << endl;

  return 0;
}

double cs(double x)
{
  double a,b,c,d;
  x = log10(x);
  a = -2.84582683378747;
  b = -1.98187129924843;
  c = -0.0528598321128271;
  d = 0.0424806763283969;
  double temp = a+(b+(c+d*x)*x)*x;
  return pow(10,temp);
}

double is(double x)
{
  double a,b,c,d;
  x = log10(x);
  double temp = a+(b+(c+d*x)*x)*x;
  return pow(10,temp);
}

double pe(double x)
{
  double a,b,c,d,e,f;
  x = log10(x);
  a = -2.35519828168804;
  b = -2.06429208226674;
  c = 0.698305851759179;
  d = 0.00502304809343589;
  e = -0.154102421097777;
  f = -0.00918861370146097;
  double temp = a+(b+(c+(d+(e+f*x)*x)*x)*x)*x;
  return pow(10,temp);
}

double npp(double x)
{
  double a,b,c,d,e,f;
  x = log10(x);
  a = -4.98224267460286;
  b = 17.0346577242561;
  c = -45.7527308923475;
  d = 68.7855291590646;
  e = -52.4572355587654;
  f = 15.8136280377002;
  double temp = a+(b+(c+(d+(e+f*x)*x)*x)*x)*x;
  return pow(10,temp);
}

double epp(double x)
{
  double a,b,c,d,e,f,g,h,i,j,k;
  x = log10(x);
  a = -5576.265643212;
  b = 90425.2750860535;
  c = -652856.065045285;
  d = 2758516.73545972;
  e = -7551838.66072215;
  f = 13997251.3011271;
  g = -17792280.3704884;
  h = 15319944.7241945;
  i = -8554643.11593555;
  j = 2798482.21068001;
  k = -407429.30716225;
  double temp = a+(b+(c+(d+(e+(f+(g+(h+(i+(j+k*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
  return pow(10,temp);
}

double mu_i(double x, double coeffs[], int nterms)
{
  x = log10(x);
  double temp=coeffs[nterms-1];
  for(int i = nterms - 2; i >=0; i--) {
    temp = coeffs[i]+temp*x;
  }
  return pow(10,temp);
}

double mu_total(double x)
{
  return mu_i(x,cs_arr,4)+mu_i(x,is_arr,4)+mu_i(x,pe_arr,6)+mu_i(x,npp_arr,6)+mu_i(x,epp_arr,11);
}
