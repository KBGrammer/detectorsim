#include "rng.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include "vector3.h"
#include <string>
#include <sstream>
#include <sys/time.h>
#include <stdlib.h>
#include "mpi.h"
#include <unistd.h>

// Last edited June 7, 2011
/********************************************************
Benchmark tests April 30, 2011, 12:32am
Times are averages, or low-hi ranges.
1 core
Events        Newton           Darth Revan      Darth Malak
1M            2.805            8.775            12.004
10M           26.137           84.664           122.34
100M          259.17                          1207

Multi core
Events        Newton(16)       Darth Revan(2)   Darth Malak(2)
1M                             9.005            18.905
10M           23.753-29.769    89.04avg         187.03
100M          236.1-244.26     877.7
500M N2       1433.5-1436.9
500M N128     1287.1-1313.4
1B N128       2215.3-2255.2
2B N128       4430.8-6667.8 (all but 2 by 4500)

Float instead of Double (old version)
Events        Newton         Darth Revan      Darth Malak
1M            5

Float instead of Double (new version) no improvement
Events        Newton         Darth Revan      Darth Malak
1M            2.689
10M           25.319
100M          

-O3 compiler option:
Events        Newton         Darth Revan      Darth Malak
10M                          32.965avg        56.777

 *******************************************************/

using namespace std;

#define huzzah endl

const double pi = acos(-1.0);
const double pi2 = 2 * pi;
const double sqrt2 = sqrt(2);
const double sqrt2o2 = sqrt2/2;
const int maxfaces = 20;
const int NDETECT = 48;
const int prec = 10;

const double r = 9.0*2.54; //22.86cm, detector radius
const double l = 15.2; // 15.2cm, detector cube size
const double slabsize = 9*2.54;
Vector3 detectors[48][8]; // Array of detector cube points
bool face_hits[48][6]; // Array of "was this face hit?"
double energydeposit[48]; // Energy deposited counter
double cosines[48]; // Cosine theta sum (then average later)
double sines[48]; // Sine theta sum (then average later)
unsigned long long first_hit[48]; // Counter of first interactions
double zm1,zm2;
double maxes[48][3]; // Array of maximum (x,y,z) values for detectors
double mins[48][3]; // Array of minimum (x,y,z) values for detectors
unsigned long long counter[48]; // Hits counter

double cs_arr[4] = {-2.84582683378747,-1.98187129924843,-0.0528598321128271,0.0424806763283969};
double is_arr[4] = {-1.27668885411202,-0.491939420500784,-0.150907754455624,0.031031257379362};
double pe_arr[6] = {-2.35519828168804,-2.06429208226674,0.698305851759179,0.00502304809343589,
		    -0.154102421097777,-0.00918861370146097};
double npp_arr[6] = {-4.98224267460286,17.0346577242561,-45.7527308923475,68.7855291590646,
		     -52.4572355587654,15.8136280377002};
double epp_arr[11] = {-5576.265643212,90425.2750860535,-652856.065045285,2758516.73545972,
		      -7551838.66072215,13997251.3011271,-17792280.3704884,15319944.7241945,
		      -8554643.11593555,2798482.21068001,-407429.30716225};

double al_sigma, al_density, al_mass, n_avogadro, lambda;
int num_discs;
double al_width, disc_width, keV;
unsigned long long numtimes;
int type;
int nprocs;
Vector3 primary_origin;

void unit_vector(Rng & random, GRng & gauss, Vector3 & vect, Vector3 & origin, int type);
void test_sources(Rng & random, GRng & gauss);
void solid_angle(Rng & random, GRng & gauss, int run, int append, int addon);
void sort_distances(double dists[], int hits[]);
void make_detector();
double dot_prod(const Vector3 & a, const Vector3 & b);
double mu_i(double x, double array[], int nterms);
double mu_total(double x);
void reset_faces();

int main (int argc, char *argv[])
{
  int rank = 0;
  int addon = 0;
  double ox,oy,oz;
  MPI::Init(argc,argv);
  rank = MPI::COMM_WORLD.Get_rank();
  nprocs = MPI::COMM_WORLD.Get_size();
  int rc;

  int t;
  if(rank == 0) {
    t = time(NULL);
  }

  rc = MPI_Bcast(&t, 1, MPI::INT, 0, MPI_COMM_WORLD);
  /* Broadcast the initial time to all other processes */

  if (rc != MPI_SUCCESS) {
    fprintf(stderr, "Oops! An error occurred in MPI_Bcast()\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //cout << t << huzzah;
  Rng random(t+rank);
  GRng gauss(t+nprocs+rank+5);
  test_sources(random,gauss);

  int scan;
  ifstream fin;
  fin.open("config.txt");
  if(argc < 8) {
    ifstream fin;
    fin.open("config.txt");
    fin >> keV >> type >> numtimes >> ox >> oy >> oz >> scan >> addon;
  }
  else {
    keV = atof(argv[1]);
    type = atoi(argv[2]);
    numtimes = atoll(argv[3]);
    ox = atof(argv[4]);
    oy = atof(argv[5]);
    oz = atof(argv[6]);
    scan = atof(argv[7]);
    addon = atoi(argv[8]);
  }

  cout << "Running " << numtimes << " at (" << ox << ',' << oy
       << ',' << oz << ") type " << type << "." << huzzah;

  primary_origin.set_coords(ox,oy,oz);

  al_sigma = 0.7*pow(10.0,-28); // cm^2
  al_density = 2.70; //  g/cm^3
  al_mass = 26.98153; // g/mol
  n_avogadro = 6.022141*pow(10.0,23);
  // Neutron mean free path
  lambda = pow(al_density/al_mass*n_avogadro*al_sigma*10000,-1);
  disc_width = 0.2; // width in cm
  num_discs = 35;
  al_width = num_discs * disc_width;

  if(scan == 1) {
    double step = 12.7*2/(sqrt(nprocs)-1);
    double step2 = step/5;
    double x0,y0;
    x0 = -12.7+rank%int(sqrt(nprocs))*step;
    y0 = int(rank/sqrt(nprocs))*step-12.7;
    int i = 0;
    int j = 0;
    ox = x0 + i * step2;
    oy = y0 + j * step2;
    oz = 22.5;
    primary_origin.set_coords(ox,oy,oz);
    solid_angle(random, gauss, rank, 1, addon);
  }
  else if(scan == 2) {
    for(int i = 0; i < 5; i++) {
      for(int j = 0; j < 5; j++) {
	oy = -12.7 + j * 2.54;
	ox = -12.7 + i * 2.54;
	oz = 22.5;
	primary_origin.set_coords(ox,oy,oz);
	solid_angle(random, gauss, rank, 1, addon);
      }
    }
  }
  else if(scan == 3) {
    ox = -12.7 + rank % 11 * 2.54;
    oy = -12.7 + int(rank/11)*2.54;
    oz = 22.5;
    primary_origin.set_coords(ox,oy,oz);
    solid_angle(random, gauss, rank, 1, addon);
  }
  else if(scan == 0) {
    solid_angle(random, gauss, rank, 0, addon);
  }

  MPI::Finalize();
  return 0;
}

void solid_angle(Rng & random, GRng & gauss, int run, int append, int addon)
{
  struct timeval tv;
  struct timezone tzone;
  gettimeofday(&tv,&tzone);
  double t0 = tv.tv_sec + tv.tv_usec/1000000.0;

  make_detector();

  double average_time = 0.0;
  ofstream fout,fout2;
  //  double x1,y1,z1;
  double energytemp0,energytemp1,dist;
  // Top of detector point, bottom point, and normal to plane
  Vector3 temptop,tempbot,norm,temp,origin,tvectop,tvectbot,temp2;
  Vector3 facea,faceb,facec,faced,facee,facef;
  facea.set_coords(0,1,0);
  faceb.set_coords(-1,0,0);
  facec.set_coords(0,0,1);
  faced.set_coords(1,0,0);
  facee.set_coords(0,0,-1);
  facef.set_coords(0,-1,0);

  double tx,ty,tz,ox,oy,oz; // photon vector and photon origin
  double x, y, z;
  bool notfound = true; // has photon hit a detector?
  int tcount; // Counter for how many detectors the photon goes through
  unsigned long long misseddetector = 0ULL; // Miss counter

  double rlength,mu; // Radiation length in cm^-1, used in keV*Exp(-dx*rlength)
  double rho = 4.51; // g/cm^3

  mu = mu_total(keV/1000);

  rlength = mu*rho;

  bool debug = false; // debug flag

  //double t0 = time(NULL);
  int mult;
  mult = 0;

  fout.open("hits.txt");
  fout2.open("misses.txt");

  double distance[maxfaces];
  int hitdetectors[maxfaces];

  fout.precision(prec);
  cout.precision(prec);
  for(unsigned int n = 0; n < numtimes; n++) {
    notfound = true;
    reset_faces();
    unit_vector(random,gauss,temp,origin,type);
    tx = temp.x; ty = temp.y; tz = temp.z;
    ox = origin.x; oy = origin.y; oz = origin.z;
    //    double tdota, tdotb, tdotc, tdotd, tdote, tdotf; // Reciprocal of temp dot faceX
    // tdota = 1/dot_prod(temp, facea); tdotb = 1/dot_prod(temp, faceb);
    // tdotc = 1/dot_prod(temp, facec); tdotd = 1/dot_prod(temp, faced);
    // tdote = 1/dot_prod(temp, facee); tdotf = 1/dot_prod(temp, facef);
    // Replace the dot products above for efficiency, since the dot products are simple
    const double tdota = 1/(temp.y);
    const double tdotb = -1/(temp.x);
    const double tdotc = 1/(temp.z);
    const double tdotd = 1/(temp.x);
    const double tdote = -1/(temp.z);
    const double tdotf = -1/(temp.y);

    int m = 0;

    for(int k = 0; k < maxfaces; k++ ) {
      distance[k] = 0.0;
      hitdetectors[k] = -1;
    }

    //if(debug) cout << tx << ' ' << ty << ' ' << tz << huzzah << huzzah;
    tcount = 0;

    while(m < 48) {
      tvectop.x = maxes[m][0]-ox; tvectop.y = maxes[m][1]-oy; tvectop.z = maxes[m][2]-l-oz;
      tvectbot.x = maxes[m][0]-l-ox; tvectbot.y = maxes[m][1]-l-oy; tvectbot.z = maxes[m][2]-oz;

      // Face a
      //dist = (tvectop.dot(facea))/temp.dot(facea);
      //dist = dot_prod(tvectop, facea) / dot_prod(temp, facea);
      if(!face_hits[m][0]) {
	//dist = dot_prod(tvectop, facea) * tdota;
	dist = tvectop.y * tdota;
	x = ox + tx*dist; y = oy + ty*dist; z = oz + tz*dist;
	//if(debug2) cout << "distance" << dist << ' ' << x << ' ' << y << ' ' << z << huzzah;
	if(x <= maxes[m][0] && x > mins[m][0] && z <= maxes[m][2] && z > mins [m][2] && dist > 0) {
	  notfound = false;
	  face_hits[m][0] = true;
	  counter[m]++;
	  distance[tcount] = dist;
	  hitdetectors[tcount] = m;
	  tcount++;
	  //if(debug) cout << m << " hit face a ";
	}
      }
      // Face b
      //dist = (tvectbot.dot(faceb))/temp.dot(faceb);
      //dist = dot_prod(tvectbot, faceb) / dot_prod(temp, faceb);
      if(!face_hits[m][1]) {
	//dist = dot_prod(tvectbot, faceb) * tdotb;
	dist = -tvectbot.x * tdotb;
	x = ox + tx*dist;y = oy + ty*dist; z = oz + tz*dist;
	//cout << "distance" << dist << ' ' << x << ' ' << y << ' ' << z << ' ' << tcount << ' ' << m << huzzah;
	if(y <= maxes[m][1] && y > mins[m][1] && z <= maxes[m][2] && z > mins [m][2] && dist > 0) {
	  //cout << counter[m] << ' ' << distance[tcount-1] << ' ' << hitdetectors[tcount-1] << huzzah;
	  notfound = false;
	  face_hits[m][1] = true;
	  counter[m]++;
	  distance[tcount] = dist;
	  hitdetectors[tcount] = m;
	  tcount++;
	  // int next = m+11;
	  // if((m==0 || m==12 || m ==24 || m==36) && !face_hits[next][3]) {
	  //   //cout << n << " Hit " << m << " face b and " << next << " face d." << endl;
	  //   face_hits[next][3] = true;
	  //   distance[tcount] = dist;
	  //   hitdetectors[tcount] = next;
	  //   tcount++;
	  //   counter[next]++;
	  // }
	  //if(debug) cout << m << " hit face b ";
	}
      }
      else {
	//cout << n << " already hit " << m << " face b." << huzzah;
      }
      // Face c
      //dist = (tvectbot.dot(facec))/temp.dot(facec);
      //dist = dot_prod(tvectbot, facec) / dot_prod(temp, facec);
      if(!face_hits[m][2]) {
	//dist = dot_prod(tvectbot, facec) * tdotc;
	dist = tvectbot.z * tdotc;
	x = ox + tx*dist;y = oy + ty*dist; z = oz + tz*dist;
	//if(debug2) cout << "distance" << dist << ' ' << x << ' ' << y << ' ' << z << huzzah;
	if(x <= maxes[m][0] && x > mins[m][0] && y <= maxes[m][1] && y > mins[m][1] && dist > 0) {
	  notfound = false;
	  face_hits[m][2] = true;
	  counter[m]++;
	  distance[tcount] = dist;
	  hitdetectors[tcount] = m;
	  tcount++;
	  // int next = m+12;
	  // if(next < 48 && !face_hits[next][4]) {
	  //   //cout << n << " Hit " << m << " face c and " << next << " face e." << endl;
	  //   face_hits[next][4] = true;
	  //   distance[tcount] = dist;
	  //   hitdetectors[tcount] = next;
	  //   tcount++;
	  //   counter[next]++;
	  // }
        //if(debug) cout << m << " hit face c ";
	}
      }
      else {
	//cout << n << " already hit " << m << " face c." << huzzah;
      }
      // Face d
      // dist = (tvectop.dot(faced))/temp.dot(faced);
      //dist = dot_prod(tvectop, faced) / dot_prod(temp, faced);
      if(!face_hits[m][3]) {
	//dist = dot_prod(tvectop, faced) * tdotd;
	dist = tvectop.x * tdotd;
	x = ox + tx*dist; y = oy + ty*dist; z = oz + tz*dist;
	//if(debug2) cout << "distance" << dist << ' ' << x << ' ' << y << ' ' << z << huzzah;
	if(y <= maxes[m][1] && y > mins[m][1] && z <= maxes[m][2] && z > mins [m][2] && dist > 0) {
	  notfound = false;
	  face_hits[m][3] = true;
	  counter[m]++;
	  distance[tcount] = dist;
	  hitdetectors[tcount] = m;
	  tcount++;
	  // int next = m-11;
	  // if((m==11 || m==23 || m ==35 || m==47) && !face_hits[next][1]) {
	  //   //cout << n << " Hit " << m << " face d and " << next << " face b." << endl;
	  //   face_hits[next][1] = true;
	  //   distance[tcount] = dist;
	  //   hitdetectors[tcount] = next;
	  //   tcount++;
	  //   counter[next]++;
	  // }
	  //if(debug) cout << m << " hit face d ";
	}
      }
      else {
	//cout << n << " already hit " << m << " face d." << huzzah;
      }
      // Face e
      //dist = (tvectop.dot(facee))/temp.dot(facee);
      //dist = dot_prod(tvectop, facee) / dot_prod(temp, facee);
      if(!face_hits[m][4]) {
	//dist = dot_prod(tvectop, facee) * tdote;
	dist = -tvectop.z * tdote;
	x = ox + tx*dist; y = oy + ty*dist; z = oz + tz*dist;
	//if(debug2) cout << "distance" << dist << ' ' << x << ' ' << y << ' ' << z << huzzah;
	if(x <= maxes[m][0] && x > mins[m][0] && y <= maxes[m][1] && y > mins [m][1] && dist > 0) {
	  notfound = false;
	  face_hits[m][4] = true;
	  counter[m]++;
	  distance[tcount] = dist;
	  hitdetectors[tcount] = m;
	  tcount++;
	  // int next = m-12;
	  // if(next > 0 && !face_hits[next][2]) {
	  //   //cout << n << " Hit " << m << " face e and " << next << " face c." << endl;
	  //   face_hits[next][2] = true;
	  //   distance[tcount] = dist;
	  //   hitdetectors[tcount] = next;
	  //   tcount++;
	  //   counter[next]++;
	  // }
	  //if(debug) cout << m << " hit face e ";
	}
      }
      else {
	//cout << n << " already hit " << m << " face e." << huzzah;
      }
      // Face f
      //dist = (tvectbot.dot(facef))/temp.dot(facef);
      //dist = dot_prod(tvectbot, facef) / dot_prod(temp, facef);
      if(!face_hits[m][5]) {
	//dist = dot_prod(tvectbot, facef) * tdotf;
	dist = -tvectbot.y * tdotf;
	x = ox + tx*dist;y = oy + ty*dist; z = oz + tz*dist;
	//if(debug2) cout << "distance" << dist << ' ' << x << ' ' << y << ' ' << z << huzzah;
	if(x <= maxes[m][0] && x > mins[m][0] && z <= maxes[m][2] && z > mins[m][2] && dist > 0) {
	  notfound = false;
	  face_hits[m][5] = true;
	  counter[m]++;
	  distance[tcount] = dist;
	  hitdetectors[tcount] = m;
	  tcount++;
	  //if(debug) cout << m << " hit face f ";
	}
      }
      else {
	//cout << n << " already hit " << m << " face f." << huzzah;
      }
      m++;
    }

    if(!notfound) { // holy shit double negative!
      sort_distances(distance, hitdetectors);
      energytemp0 = keV;
      //if(debug) cout << tcount << huzzah;
      first_hit[hitdetectors[tcount-1]]++;
      //cosine[hitdetectors[p]] += 
      //cout << tcount << '\t';
      for(int p = maxfaces - 1; p > 0; p-=2) {
	if(hitdetectors[p] != -1) {
	  //cout << p << '\t';
	  if(p == tcount - 1) {
	    //cout << "yes";
	    //double c = (oy+ty*distance[p])/distance[p];
	    //double s = (ox+tx*distance[p])/distance[p];
	    // cout << ox << '\t' << oy << '\t' << oz << '\t' << distance[p] << '\t';
	    // cout << tx << '\t' << ty << '\t' << tz << '\t'; 
	    cosines[hitdetectors[tcount-1]] += ty;
	    sines[hitdetectors[tcount-1]] += tx;
	    //cout << c << '\t' << s << huzzah;
	    // double cossin = (c*c+s*s);
	    // cout << c << '\t' << s << '\t' << cossin << huzzah;
	  }
	  energytemp1 = energytemp0*exp((distance[p]-distance[p-1])*rlength);
	  energydeposit[hitdetectors[p]] += energytemp0-energytemp1;
	  energytemp0 = energytemp1;
	}
      }
      //cout << huzzah;
    }
    if(notfound) {
      misseddetector++;
    }

    // Prints a vector to a file every once in a while, since gnuplot can be overwhelmed with points.
    if(n < 10000 && run == 0) {
      if(notfound) {
	//temp.return_coords(tx,ty,tz);
	if(tz > 0) {
	  norm.set_coords(0,0,1);
	  temp2.set_coords(0-ox,0-oy,2*l-oz);
	}
	else {
	  norm.set_coords(0,0,-1);
	  temp2.set_coords(0-ox,0-oy,-2*l-oz);
	}
	// get distance to intersection
	dist = dot_prod(temp2, norm)/dot_prod(temp, norm);
	x = ox + temp.x*dist;
	y = oy + temp.y*dist;
	z = oz + temp.z*dist;
	fout2 << x << ' ' << y << ' ' << z << huzzah;
      }
      else {
	if (debug) cout << huzzah;
	for(int p = maxfaces - 1; p > 0; p-=2) {
	  if(hitdetectors[p] != -1) {
	    double d0 = distance[p];
	    fout << hitdetectors[tcount-1] << ' ' << ox + temp.x*d0 << ' '
		 << oy + temp.y*d0 << ' ' << oz + temp.z*d0 << huzzah;
	    d0 = distance[p-1];
	    fout << hitdetectors[tcount-1] << ' ' << ox + temp.x*d0 << ' '
		 << oy + temp.y*d0 << ' ' << oz + temp.z*d0 << huzzah;
	  }
	}
      }
    }
  }
    
  // Finally, output the result counts.
  fout.close();
  
  string filename;
  ostringstream oss;
  oss << "counts" << nprocs << "x" << keV << "-"
      << type << "-" << double(numtimes) << "avg.txt";
  filename += oss.str();

  if(append == 1 || addon == 1) {
    fout.open(filename.c_str(), ios_base::app);
  }
  else {
    fout.open(filename.c_str());
  }
    
  ofstream fout3;
  string filename2;
  ostringstream oss2;
  oss2 << "counts" << nprocs << "x" << keV << "-"
       << type << "-" << double(numtimes) << ".txt";
  filename2 += oss2.str();
    
  double pox,poy,poz;
  primary_origin.return_coords(pox,poy,poz);
  double orig[3]={pox,poy,poz};
  gettimeofday(&tv,&tzone);
  double tf = tv.tv_sec + tv.tv_usec/1000000.0 - t0;
  //double tf = time(NULL) - t0;
  cout.precision(5);
  cout << "Run " << run << " took " << tf << " seconds." << huzzah;
  cout.precision(prec);
  
  if(run != 0) { // Send data to rank 0 processor
    MPI_Send(&energydeposit, NDETECT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&cosines, NDETECT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&sines, NDETECT, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&counter, NDETECT, MPI::UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&first_hit, NDETECT, MPI::UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&misseddetector, 1, MPI::UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&orig, 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    MPI_Send(&tf, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }
  else {
    if(append == 1 || addon == 1) {
      fout3.open(filename2.c_str(), ios_base::app);
    }
    else {
      fout3.open(filename2.c_str());
    }
    fout3.precision(prec);
    fout.precision(prec);
    double torig[3];
    double time = 0;
    for(int i = 0; i < nprocs; i++) {
      MPI_Status status;
      double ted[NDETECT];
      double tcos[NDETECT];
      double tsin[NDETECT];
      unsigned long long tfh[NDETECT];
      unsigned long long tco[NDETECT];
      unsigned long long tmis;
      if(i > 0) { // Receive data from other processors
	MPI_Recv(&ted, NDETECT, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&tcos, NDETECT, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&tsin, NDETECT, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&tco, NDETECT, MPI::UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&tfh, NDETECT, MPI::UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&tmis, 1, MPI::UNSIGNED_LONG, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&torig, 3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
	MPI_Recv(&time, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);

      }
      else if(i == 0) {
	for(int j = 0; j < NDETECT; j++) {
	  ted[j] = energydeposit[j];
	  tco[j] = counter[j];
	  tfh[j] = first_hit[j];
	}
	tmis = misseddetector;
	torig[0] = pox;
	torig[1] = poy;
	torig[2] = poz;
      }
      average_time += time;
      fout3 << "R " << i << " This is run " << i << huzzah;
      fout3 << "R " << i << " Detector Num Total FirstHits "
	    << "FirstFraction TotalEnergy EnergyFraction x y z cos() sin()" << huzzah;
      unsigned long long sum1 = 0;
      double sum2 = 0;
      for(int j = 0; j < NDETECT; j++) {
	sum1 += tfh[j];
	sum2 += ted[j];
	fout3 << "R " << i << ' ' << "Detector " << j << ' ' << tco[j]/2 << ' ';
	fout3 << tfh[j] << ' ' << tfh[j]/double(numtimes);
	fout3 << ' ' << ted[j] << ' ';
	fout3 << ted[j]/(double(numtimes)*keV);
	fout3 << ' ' << torig[0] << ' ' << torig[1] << ' ' << torig[2] << ' '
	      << tcos[j]/tfh[j] << ' ' << tsin[j]/tfh[j] << huzzah;
	if(i!= 0) {
	  first_hit[j] += tfh[j];
	  energydeposit[j] += ted[j];
	  cosines[j] += tcos[j];
	  sines[j] += tsin[j];
	  counter[j] += tco[j];
	}
      }
      fout3 << "R " << i << ' ' << "Detector Hits " << sum1 << ' ';
      fout3 << tmis << ' ' << tmis+sum1 << huzzah;
      fout3 << "R " << i << ' ' << "Energy Deposited " << sum2;
      fout3 << " Total Energy " <<  double(keV*numtimes) << huzzah;
      fout3 << "R " << i << ' ' << "Total Fraction of Energy Deposited: ";
      fout3 << sum2 / (keV*numtimes) << " = ";
      fout3 << sum2 / (keV*numtimes)*4 << "*pi." << huzzah;
      fout3 << "R " << i << ' ' << "Total Fraction of First Hits: ";
      fout3 << double(sum1) / double(numtimes) << " = ";
      fout3 << double(sum1) / double(numtimes)*4 << "*pi." << huzzah;
      fout3 << "R " << i << " The photon energy was " << keV;
      fout3 << " keV and there were " << numtimes << " events." << huzzah;
      fout3 << "R " << i << " The run type was " << type << " and the origin was (";
      fout3 << torig[0] << ' ' << torig[1] << ' ' << torig[2] << ")." << huzzah;
      //cout << torig[0] << ' ' << torig[1] << ' ' << torig[2] << huzzah;
      if(i != 0) {
	misseddetector += tmis;
      }
    }
    if(run == 0) {
      fout << "Detector Num Total FirstHits "
	   << "FirstFraction TotalEnergy EnergyFraction x y z cos() sin()" << huzzah;
      unsigned long long sum3 = 0;
      double sum4 = 0;
      for(int e = 0; e < NDETECT; e++) {
	sum3 += first_hit[e];
	sum4 +=  energydeposit[e];
	fout << "Detector " << e << ' ' << counter[e]/2 << ' ';
	fout << first_hit[e] << ' ' << first_hit[e]/double(numtimes*nprocs);
	fout << ' ' << energydeposit[e] << ' ';
	fout << energydeposit[e]/(double(numtimes)*keV*nprocs);
	fout << ' ' << torig[0] << ' ' << torig[1] << ' ' << torig[2] << ' ' 
	     << cosines[e]/first_hit[e] << ' ' << sines[e]/first_hit[e] << huzzah;
      }
      fout << "Detector Hits " << sum3 << ' ' << misseddetector;
      fout << ' ' << misseddetector+sum3 << huzzah;
      fout << "Energy Deposited " << sum4 << " Total Energy ";
      fout << double(keV*numtimes*nprocs) << huzzah;
      fout << "Total Fraction of Energy Deposited: ";
      fout << sum4 / double(keV*numtimes*nprocs) << " = ";
      fout << sum4 / double(keV*numtimes*nprocs)*4 << "*pi." << huzzah;
      fout << "Total Fraction of First Hits: ";
      fout << double(sum3) / double(numtimes*nprocs) << " = ";
      fout << double(sum3) / double(numtimes*nprocs)*4 << "*pi." << huzzah;
      fout << "The photon energy was " << keV << " keV and there were ";
      fout << numtimes*nprocs << " events." << huzzah;
      fout << "The run type was " << type << " and the origin was (";
      fout << torig[0] << ' ' << torig[1] << ' ' << torig[2] << ")." << huzzah;
    }
    cout << "Done." << huzzah;
    average_time = (average_time + tf) / nprocs;
    gettimeofday(&tv,&tzone);
    tf = tv.tv_sec + tv.tv_usec/1000000.0 - t0;
    cout.precision(5);
    fout.precision(5);
    cout << "Operation took " << tf << " seconds (average " 
	 << average_time << ") on " << nprocs << " cores." << huzzah;
    fout << "Operation took " << tf << " seconds (average " 
	 << average_time << ") on " << nprocs << " cores." << huzzah;
    fout.close();
    fout2.close();
    fout3.close();
  }
}

void unit_vector(Rng & random, GRng & gauss, Vector3 &vect, Vector3 & origin, int type)
{
  double u,theta,temp,x,y,z,ox,oy,oz;
  u = random.rdouble()*2.0-1.0;
  theta = random.rdouble_exc()*pi2;
  temp = sqrt(1-u*u);
  x = temp*cos(theta);
  y = temp*sin(theta);
  z = u;
  //vect.set_coords(x,y,z);
  vect.x = x;
  vect.y = y;
  vect.z = z;
  // double mag1 = vect.magnitude();
  // vect.make_unitvector();
  // double mag2 = vect.magnitude();
  // cout << mag2 - mag1 << huzzah;

  //  double slabsize = 9*2.54;

  ox = primary_origin.x;
  oy = primary_origin.y;
  oz = primary_origin.z;

  if(type == 1) { // Point source at origin
    origin.set_coords(ox,oy,oz);
  }
  else if(type == 2) {// cylinder source
    double tx,ty;
    do{
      tx = (random.rdouble()*2.0-1.0);
      ty = (random.rdouble()*2.0-1.0);
    }while(tx*tx+ty*ty>1);
    tx = tx*5.715; // 4.5*2.54/2 = radius in cm
    ty = ty*5.715;
    double tz = (random.rdouble()-0.5)*0.22; // width
    origin.set_coords(ox+tx,oy+ty,oz+tz);
  }
  else if(type == 3) { // pill source
    origin.set_coords(ox,oy,oz+double(random.rdouble()*60.0-30.0));
  }
  else if(type == 4) { // 1x1x1cm cube at origin
    origin.set_coords(ox+random.rdouble()-0.5,oy+random.rdouble()-0.5,oz+random.rdouble()-0.5);
  }
  else if(type == 5) { // 10x10x1cm slab at origin
    origin.set_coords(ox+random.rdouble()*slabsize-slabsize/2.0,oy+random.rdouble()*slabsize-slabsize/2.0,random.rdouble()-0.5+oz);
  }
  else if(type == 6) { // 10x10x1cm slab at origin rotated +pi/4
    x = random.rdouble()*slabsize-slabsize/2.0;
    y = random.rdouble()*slabsize-slabsize/2.0;
    z = random.rdouble()-0.5;
    origin.set_coords(x+ox,oy+sqrt2o2*(y-z),oz+sqrt2o2*(y+z));
  }
  else if(type == 7) { // 10x10x1cm slab at origin rotated +pi/4
    x = random.rdouble()*slabsize-slabsize/2.0;
    y = random.rdouble()*slabsize-slabsize/2.0;
    z = random.rdouble()-0.5;
    origin.set_coords(ox+sqrt2o2*(x+z),y+oy,oz+sqrt2o2*(-x+z));
  }
  else if(type == 8) { // 10x10x0cm sheet at origin
    x = random.rdouble()*10.0-5.0;
    y = random.rdouble()*10.0-5.0;
    z = 0;
    origin.set_coords(x+ox,y+oy,z+oz);
  }
  else if(type == 9) { // 10x10x0cm sheet at origin rotated about y-axis +pi/4
    x = random.rdouble()*10.0-5.0;
    y = random.rdouble()*10.0-5.0;
    z = random.rdouble()-0.5;
    origin.set_coords(sqrt2o2*(x-z)+ox,y+oy,sqrt2o2*(x+z)+oz);
  }
  else if(type == 10) {
    do {
      z = -lambda*log(random.rdouble());
    } while(z > al_width);
    int disc = z/disc_width;
    z = z - disc * disc_width;
    z = z + disc * disc_width * 2;
    double r = 4.0*sqrt(random.rdouble());
    double theta = pi2*random.rdouble();
    x = r * cos(theta);
    y = r * sin(theta);
    //    z = random.rdouble()*10.0-5.0;
    origin.set_coords(ox+x,oy+y,oz+z);
  }
}

// insertion sort algorithm
// from http://mathbits.com/mathbits/compsci/arrays/Insertion.htm
void sort_distances(double dists[], int hits[])
{
  int i, j;
  double key;
  char ckey;
  for(j = 1; j < maxfaces && dists[j] != -1; j++) {  // Start with 1 (not 0)
    key = fabs(dists[j]);
    ckey = hits[j];
    // Smaller values move up
    for(i = j - 1; (i >= 0) && (fabs(dists[i]) < key); i--) {
      dists[i+1] = dists[i];
      hits[i+1] = hits[i];
    }
    dists[i+1] = key;    //Put key into its proper location
    hits[i+1] = ckey;
  }
  for(int l = 0; l < maxfaces && hits[l-1] != -1; l++) {
    if(l%2 == 0 && hits[l] != hits[l+1]) {
      hits[l+2] = hits[l+1];
      hits[l+1] = hits[l];
    }
  }
  return;
}

void test_sources(Rng & random, GRng & gauss)
{
  ofstream fout;
  Vector3 temp,origin;
  double x,y,z,tx,ty,tz,ox,oy,oz;
  tx=ty=tz=0.0;
  for(int i = 12; i < 11; i++ ) {
    if(i == 1) {
      fout.open("plotfiles/pointsource.txt");
    }
    if(i == 2) {
      fout.open("plotfiles/linesource.txt");
    }
    if(i == 3) {
      fout.open("plotfiles/pillsource.txt");
    }
    if(i == 4) {
      fout.open("plotfiles/cubesource.txt");
    }
    if(i == 5) {
      fout.open("plotfiles/slabsource.txt");
    }
    if(i == 6) {
      fout.open("plotfiles/slabrotsourcex.txt");
    }
    if(i == 7) {
      fout.open("plotfiles/slabrotsourcey.txt");
    }
    if(i == 8) {
      fout.open("plotfiles/sheetsource.txt");
    }
    if(i == 9) {
      fout.open("plotfiles/sheetrotsourcey.txt");
    }
    if(i == 10) {
      fout.open("plotfiles/al_source.txt");
    }
    for(int j = 0; j < 20000; j++) {
      unit_vector(random,gauss,temp,origin,i);
      temp.return_coords(tx,ty,tz);
      origin.return_coords(ox,oy,oz);
      //x = tx+ox;y=ty+oy;z=tz+oz;
      x=ox;y=oy;z=oz;
      fout << x << ' ' << y << ' ' << z << huzzah;
    }
    fout.close();
  }
}

void make_detector()
{
  ofstream fout;
  fout.open("array.txt");
  // This for loop creates the detector array as a matrix of vectors.
  // It also output coordinates to a file for plotting.
  for(int j = 0; j < 48; j++) {
    energydeposit[j] = 0.0;
    first_hit[j] = 0;
    double i = int(j/12);
    zm2 = (i-2)*l;
    zm1 = zm2+l;
    maxes[j][2] = zm1;
    

    // Create the detector array
    switch (j%12) {
      //0, top right inside, clockwise, 4 top right, outside, clockwise
    case 0:
      detectors[j][0].set_coords(l,r,zm2);
      detectors[j][1].set_coords(l,r,zm1);
      detectors[j][2].set_coords(0.0,r,zm1);
      detectors[j][3].set_coords(0.0,r,zm2);
      detectors[j][4].set_coords(l,r+l,zm2);
      detectors[j][5].set_coords(l,r+l,zm1);
      detectors[j][6].set_coords(0.0,r+l,zm1);
      detectors[j][7].set_coords(0.0,r+l,zm2);
      maxes[j][0] = l; maxes[j][1] = r+l;
      break;
    case 1: // corner
      detectors[j][0].set_coords(l,l,zm2);
      detectors[j][1].set_coords(l,l,zm1);
      detectors[j][2].set_coords(l+l,l,zm1);
      detectors[j][3].set_coords(l+l,l,zm2);
      detectors[j][4].set_coords(l,l+l,zm2);
      detectors[j][5].set_coords(l,l+l,zm1);
      detectors[j][6].set_coords(l+l,l+l,zm1);
      detectors[j][7].set_coords(l+l,l+l,zm2);
      maxes[j][0] = l+l; maxes[j][1] = l+l;
      break;
    case 2:
      detectors[j][0].set_coords(r,0.0,zm2);
      detectors[j][1].set_coords(r,0.0,zm1);
      detectors[j][2].set_coords(r,l,zm1);
      detectors[j][3].set_coords(r,l,zm2);
      detectors[j][4].set_coords(r+l,0.0,zm2);
      detectors[j][5].set_coords(r+l,0.0,zm1);
      detectors[j][6].set_coords(r+l,l,zm1);
      detectors[j][7].set_coords(r+l,l,zm2);
      maxes[j][0] = r+l; maxes[j][1] = l;
      break;
    case 3:
      detectors[j][0].set_coords(r,-l,zm2);
      detectors[j][1].set_coords(r,-l,zm1);
      detectors[j][2].set_coords(r,0.0,zm1);
      detectors[j][3].set_coords(r,0.0,zm2);
      detectors[j][4].set_coords(r+l,-l,zm2);
      detectors[j][5].set_coords(r+l,-l,zm1);
      detectors[j][6].set_coords(r+l,0.0,zm1);
      detectors[j][7].set_coords(r+l,0.0,zm2);
      maxes[j][0] = r+l; maxes[j][1] = 0;
      break;
    case 4: // corner
      detectors[j][0].set_coords(l,-l,zm2);
      detectors[j][1].set_coords(l,-l,zm1);
      detectors[j][2].set_coords(l+l,-l,zm1);
      detectors[j][3].set_coords(l+l,-l,zm2);
      detectors[j][4].set_coords(l,-l-l,zm2);
      detectors[j][5].set_coords(l,-l-l,zm1);
      detectors[j][6].set_coords(l+l,-l-l,zm1);
      detectors[j][7].set_coords(l+l,-l-l,zm2);
      maxes[j][0] = l+l; maxes[j][1] = -l;
      break;
    case 5:
      detectors[j][0].set_coords(l,-r,zm2); // Top right inside
      detectors[j][1].set_coords(l,-r,zm1); // Bottom left inside
      detectors[j][2].set_coords(0.0,-r,zm1); // Bottom right inside
      detectors[j][3].set_coords(0.0,-r,zm2); // Top left inside
      detectors[j][4].set_coords(l,-r-l,zm2); // Top right outside
      detectors[j][5].set_coords(l,-r-l,zm1); // Bottom left outside
      detectors[j][6].set_coords(0.0,-r-l,zm1); // Bottom right outside
      detectors[j][7].set_coords(0.0,-r-l,zm2); // Top left outside
      maxes[j][0] = l; maxes[j][1] = -r;
      break;
    case 6:
      detectors[j][0].set_coords(0.0,-r,zm2); //Top Right, inside
      detectors[j][1].set_coords(0.0,-r,zm1);
      detectors[j][2].set_coords(-l,-r,zm1);
      detectors[j][3].set_coords(-l,-r,zm2);
      detectors[j][4].set_coords(0.0,-r-l,zm2);
      detectors[j][5].set_coords(0.0,-r-l,zm1);
      detectors[j][6].set_coords(-l,-r-l,zm1);
      detectors[j][7].set_coords(-l,-r-l,zm2);
      maxes[j][0] = 0; maxes[j][1] = -r;
      break;
    case 7: // corner
      detectors[j][0].set_coords(-l,-l,zm2);
      detectors[j][1].set_coords(-l,-l,zm1);
      detectors[j][2].set_coords(-l-l,-l,zm1);
      detectors[j][3].set_coords(-l-l,-l,zm2);
      detectors[j][4].set_coords(-l,-l-l,zm2);
      detectors[j][5].set_coords(-l,-l-l,zm1);
      detectors[j][6].set_coords(-l-l,-l-l,zm1);
      detectors[j][7].set_coords(-l-l,-l-l,zm2);
      maxes[j][0] = -l; maxes[j][1] = -l;
      break;
    case 8:
      detectors[j][0].set_coords(-r,-l,zm2);
      detectors[j][1].set_coords(-r,-l,zm1);
      detectors[j][2].set_coords(-r,0.0,zm1);
      detectors[j][3].set_coords(-r,0.0,zm2);
      detectors[j][4].set_coords(-l-r,-l,zm2);
      detectors[j][5].set_coords(-l-r,-l,zm1);
      detectors[j][6].set_coords(-l-r,0.0,zm1);
      detectors[j][7].set_coords(-l-r,0.0,zm2);
      maxes[j][0] = -r; maxes[j][1] = 0;
      break;
    case 9:
      detectors[j][0].set_coords(-r,0.0,zm2);
      detectors[j][1].set_coords(-r,0.0,zm1);
      detectors[j][2].set_coords(-r,l,zm1);
      detectors[j][3].set_coords(-r,l,zm2);
      detectors[j][4].set_coords(-l-r,0.0,zm2);
      detectors[j][5].set_coords(-l-r,0.0,zm1);
      detectors[j][6].set_coords(-l-r,l,zm1);
      detectors[j][7].set_coords(-l-r,l,zm2);
      maxes[j][0] = -r; maxes[j][1] = l;
      break;
    case 10: // corner
      detectors[j][0].set_coords(-l,l,zm2);
      detectors[j][1].set_coords(-l,l,zm1);
      detectors[j][2].set_coords(-l-l,l,zm1);
      detectors[j][3].set_coords(-l-l,l,zm2);
      detectors[j][4].set_coords(-l,l+l,zm2);
      detectors[j][5].set_coords(-l,l+l,zm1);
      detectors[j][6].set_coords(-l-l,l+l,zm1);
      detectors[j][7].set_coords(-l-l,l+l,zm2);
      maxes[j][0] = -l; maxes[j][1] = l+l;
      break;
    case 11:
      detectors[j][0].set_coords(0.0,r,zm2);
      detectors[j][1].set_coords(0.0,r,zm1);
      detectors[j][2].set_coords(-l,r,zm1);
      detectors[j][3].set_coords(-l,r,zm2);
      detectors[j][4].set_coords(0.0,r+l,zm2);
      detectors[j][5].set_coords(0.0,r+l,zm1);
      detectors[j][6].set_coords(-l,r+l,zm1);
      detectors[j][7].set_coords(-l,r+l,zm2);
      maxes[j][0] = 0; maxes[j][1] = r+l;
      break;
    }
    mins[j][0] = maxes[j][0] - l;
    mins[j][1] = maxes[j][1] - l;
    mins[j][2] = maxes[j][2] - l;
    counter[j] = 0;
    for(int d = 0; d < 8; d++) {
      fout << j << ' ' << detectors[j][d].x << ' ' << detectors[j][d].y << ' '
	   << detectors[j][d].z << huzzah;
    }
  }
  fout.close();
}

double dot_prod(const Vector3 & a, const Vector3 & b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
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
/*  if(keV == 662) mu = 0.07754;
  if(keV == 2200) mu = 0.03997;
  if(keV == 1000) mu = 0.05848;
  if(keV == 488) mu = 0.1004;
  */
  return mu_i(x,cs_arr,4)+mu_i(x,is_arr,4)+mu_i(x,pe_arr,6)+mu_i(x,npp_arr,6)+mu_i(x,epp_arr,11);
}

void reset_faces()
{
  for(int i = 0; i < 48; i++) {
    for(int j = 0; j < 6; j++) {
      face_hits[i][j] = false;
    }
  }
}
