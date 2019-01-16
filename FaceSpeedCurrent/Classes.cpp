/**
	Adds the definitions of the functions for the prototypes created in the Class header file
*/

#include "Classes.hh"

double volume::cal_flux_sum(double fluxN,double fluxS,double fluxE,double fluxW){
	return (fluxW*n_W)+(fluxE*n_E)+(fluxS*n_S)+(fluxN*n_N);
	}

	
double face::cal_dif_flux_on_face(const vector<double> &T, double alpha){
		return A*alpha*(T[Neigh_2_ID]-T[Neigh_1_ID])/2;
	}