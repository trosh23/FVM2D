#ifndef CLASSES
#define CLASSES

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

class volume{
	public:

		int ID;

		int Neigh_W_ID;
		int Neigh_E_ID;
		int Neigh_N_ID;
		int Neigh_S_ID;

		double n_W;
		double n_E;
		double n_N;
		double n_S;

		double V;
		double midpointx;
		double midpointy;
		int bound_flag;

		double Temp;

		double fluxsum;	
		double cal_flux_sum(double,double,double,double);
	};


class face{
	public:

		int ID;

		int Neigh_1_ID;
		int Neigh_2_ID;

		double n;
		double A;

		int bound_flag;
		double flux;

		double cal_dif_flux_on_face(const vector<double>&,double);
	};



#endif