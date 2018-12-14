#include "Classes.cpp"


int main(){

			cout << endl << endl;
			cout << "Enter 1 for custom parameters, or any other value for default parameters:  ";

			int inputs;
			cin >> inputs;

			int M;
			int N;
			double Lx;
			double Ly;
			double alpha;
			double source;

					switch(inputs){
						case 1:
							cout << endl << endl;
							cout << "Number of cells in x-direction:  ";
							cin >> M;

							cout << endl << endl;
							cout << "Number of cells in y-direction:  ";
							cin >> N;

							cout << endl << endl;
							cout << "Length in x-direction:  ";
							cin >> Lx;

							cout << endl << endl;
							cout << "Length in y-direction:  ";
							cin >> Ly;

							cout << endl << endl;
							cout << "Diffusivity value:  ";
							cin >> alpha;

							cout << endl << endl;
							cout << "Source term value:  ";
							cin >> source;
						
							cout << endl << endl;
						break;

						default:
							M = 10;
							N = 10;
							Lx = 1;
							Ly = 2;
							alpha = 0.1;
							source = -0.001;
						break;
					}

	double delx = Lx/M;
	double dely = Ly/N;

	vector<volume> volume_array;
	vector<face> face_array;
	vector<double> Temp;

	volume_array.reserve(M*N);
	Temp.reserve(M*N);
	face_array.reserve(2*M*N);

	int k;
	int f;

	for(int j = 0;j<N;j++){
		for(int i = 0;i<M;i++){
			k = ((i*N)+j);
			
			//Initializing values for control volumes
			volume_array[k].ID = k;

			volume_array[k].midpointx = (i*delx) + (delx/2);
			volume_array[k].midpointy = (j*dely) + (dely/2);

			volume_array[k].n_N = 1;
			volume_array[k].n_S = -1;
			volume_array[k].n_E = 1;
			volume_array[k].n_W = -1;

			volume_array[k].V = delx*dely;

			if(j == 0)			volume_array[k].bound_flag = 1;
			else if(j == N-1)	volume_array[k].bound_flag = 1;

			if(i== 0)			volume_array[k].bound_flag = 1;
			else if(i == M-1)	volume_array[k].bound_flag = 1;


		}
	}

	int boundarycondition;

	cout << endl << endl;
	cout << "Enter 1 for pre-defined boundary conditions, or anything else for the default periodic boundary conditions:  ";
	cin >> boundarycondition;

	switch(boundarycondition){
		case 1:
			//Define face IDs for non-periodic
			for(int j = 0;j<N;j++){
				for(int i = 0;i<M;i++){
					k = ((j*N)+i);
					f = k+(M*N);
					face_array[k].ID = k;
					face_array[f].ID = f;
					face_array[k].A = dely;
					face_array[f].A = delx;

					/*
					if(j == N-1) 	volume_array[k].Neigh_N_ID = -4;
					else 			volume_array[k].Neigh_N_ID = k+M;

					if(j == 0) 		volume_array[k].Neigh_S_ID = -2;
					else 			volume_array[k].Neigh_S_ID = k-M;

					if(i == M-1) 	volume_array[k].Neigh_E_ID = -3;
					else 			volume_array[k].Neigh_E_ID = k+1;

					if(i == 0) 		volume_array[k].Neigh_W_ID = -1;
					else 			volume_array[k].Neigh_W_ID = k-1;
					*/

					if(j == N-1) 	face_array[f].Neigh_2_ID = -4;
					else 			face_array[f].Neigh_2_ID = k;

					if(j == 0) 		face_array[f].Neigh_1_ID = -2;
					else 			face_array[f].Neigh_1_ID = k-M;

					if(i == M-1) 	face_array[k].Neigh_2_ID = -3;
					else 			face_array[k].Neigh_2_ID = k;

					if(i == 0) 		face_array[k].Neigh_1_ID = -1;
					else 			face_array[k].Neigh_1_ID = k-1;
				}
			}

		break;

		default:
			face_array.reserve(2*M*N);

			//Define face IDs for periodic
			for(int j = 0;j<N;j++){
				for(int i = 0;i<M;i++){
					k = ((j*N)+i);

					
					f = k+(M*N);
					face_array[k].ID = k;
					face_array[f].ID = f;
					face_array[k].A = dely;
					face_array[f].A = delx;

					if(i == 0) 		face_array[k].Neigh_1_ID = (M-1)+(M*j);
					else 			face_array[k].Neigh_1_ID = k-1;
					
					if(i == (M-1)) 	face_array[k].Neigh_2_ID = (j*M);
					else 			face_array[k].Neigh_2_ID = k;
					
					if(j == 0) 		face_array[f].Neigh_1_ID = ((N-1)*N)+i;
					else 			face_array[f].Neigh_1_ID = k-M;
					
					if(j == (N-1)) 	face_array[f].Neigh_2_ID = i;
					else 			face_array[f].Neigh_2_ID = k;
					

				/*
						if(j == N-1) 	volume_array[k].Neigh_N_ID = i;
						else 			volume_array[k].Neigh_N_ID = k+M;

						if(j == 0) 		volume_array[k].Neigh_S_ID = ((M-1)*N) + i;
						else 			volume_array[k].Neigh_S_ID = k-M;

						if(i == M-1) 	volume_array[k].Neigh_E_ID = j*N;
						else 			volume_array[k].Neigh_E_ID = k+1;

						if(i == 0) 		volume_array[k].Neigh_W_ID = (j*N)+(M-1);
						else 			volume_array[k].Neigh_W_ID = k-1;
				*/
					}
				}
		break;
	}

	double timer;

	cout << endl << endl;
	cout << "Enter runtime for the simulation in seconds:  ";
	cin >> timer;

		clock_t timers;
		timers = clock();


	double delt = 0.1*(0.25*delx*dely);

	double denom = sqrt(pow((Lx/2),2)+pow((Lx/2),2));
	double numa;
	double numb;

	ofstream Tout;
  	Tout.open ("Output.txt");
  	
	for(int j = 0;j<N;j++){
		for(int i = 0;i<M;i++){
			k = ((j*N)+i);

			numa = (volume_array[k].midpointx - (Lx/2));
			numb = (volume_array[k].midpointy - (Ly/2));

			numa = pow(numa,2);
			numb = pow(numb,2);

			Temp[k] = (numa + numb)/denom;

			Tout << Temp[k] << " ";
		}

		Tout << endl;
	}


	int logcount = 0;
	int logger = timer/(50*delt);

	double fluxN;
	double fluxS;
	double fluxE;
	double fluxW;

	double fluxsum;


	for(double t = delt;t<timer;t=t+delt){

		for(int j = 0;j<N;j++){
			for(int i = 0;i<M;i++){
				k = ((j*N)+i);
				f = k+(M*N);

				if(boundarycondition == 1){
					if(j == N-1) 	fluxN = 0;
					else 			fluxN = face_array[f+M].cal_dif_flux_on_face(Temp,alpha);

					if(j == 0) 		fluxS = ((delx*alpha)*((Temp[k]-1)/2));
					else 			fluxS = face_array[f].cal_dif_flux_on_face(Temp,alpha);

					if(i == M-1) 	fluxE = ((dely*alpha)*((-1-(Temp[k])/2)));
					else 			fluxE = face_array[k+1].cal_dif_flux_on_face(Temp,alpha);

					if(i == 0) 		fluxW = -abs(source/3);
					else 			fluxW = face_array[k].cal_dif_flux_on_face(Temp,alpha);
				}

				else{
					fluxW = face_array[k].cal_dif_flux_on_face(Temp,alpha);
					fluxE = face_array[k+1].cal_dif_flux_on_face(Temp,alpha);
					fluxS = face_array[f].cal_dif_flux_on_face(Temp,alpha);
					if(j == (N-1)) fluxN = face_array[i].cal_dif_flux_on_face(Temp,alpha);
					else 		   fluxN = face_array[f+M].cal_dif_flux_on_face(Temp,alpha);
				}

				
				volume_array[k].fluxsum = volume_array[k].cal_flux_sum(fluxN,fluxS,fluxE,fluxW);

				/*
				if(boundarycondition == 1){
					if(j == N-1) 	fluxN = 0;
					else 			fluxN = ((delx*alpha)*((volume_array[volume_array[k].Neigh_N_ID].Temp - volume_array[k].Temp)/2));


					if(j == 0) 		fluxS = ((delx*alpha)*((volume_array[k].Temp-1)/2));
					else 			fluxS = ((delx*alpha)*((volume_array[k].Temp - volume_array[volume_array[k].Neigh_S_ID].Temp)/2));

					if(i == M-1) 	fluxE = ((delx*alpha)*((-1-(volume_array[k].Temp)/2)));
					else 			fluxE = ((dely*alpha)*((volume_array[volume_array[k].Neigh_E_ID].Temp - volume_array[k].Temp)/2));

					if(i == 0) 		fluxW = -abs(source/3);
					else 			fluxW = ((dely*alpha)*((volume_array[k].Temp - volume_array[volume_array[k].Neigh_W_ID].Temp)/2));
				}

				else{
					fluxW = ((dely*alpha)*((volume_array[k].Temp - volume_array[volume_array[k].Neigh_W_ID].Temp)/2));
					fluxE = ((dely*alpha)*((volume_array[volume_array[k].Neigh_E_ID].Temp - volume_array[k].Temp)/2));
					fluxN = ((delx*alpha)*((volume_array[volume_array[k].Neigh_N_ID].Temp - volume_array[k].Temp)/2));
					fluxS = ((delx*alpha)*((volume_array[k].Temp - volume_array[volume_array[k].Neigh_S_ID].Temp)/2));
				}

				volume_array[k].fluxsum = fluxE-fluxW+fluxN-fluxS;
				*/

				Temp[k]= Temp[k] + ((delt/(delx*dely))*volume_array[k].fluxsum) + (delt*source);
			}
		}


		
		if(logcount == logger){
			Tout << endl << endl;
			Tout << "Current Time:  " << t;
			Tout << endl << "Data: " << endl << endl; 
			for(int j = N-1 ;j>=0;j--){
				for(int i = 0;i<M;i++){
					k = ((j*N)+i);
					Tout << Temp[k] << " ";
				}

				Tout << endl;
			}
			
			logcount = 0;
		}
		else logcount = logcount+1;

	}

	timers = clock() - timers;

	cout << endl << endl;
			
	cout << "\nTime elapsed: "<< (float)timers/CLOCKS_PER_SEC << " seconds" << endl << endl;

	cout << endl << endl;

}