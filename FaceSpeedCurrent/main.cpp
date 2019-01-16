/**
	Main loop which includes the definition of the inputs, face definitions and solution loop
*/

#include "Classes.cpp"	//Including source file which includes additional headers

int main(){

			//Prompt for input to determine if the user wishes to use custom parameters, or to use the hard-coded defaults

			cout << endl << endl;
			cout << "Enter 1 for custom parameters, or any other value for default parameters:  ";

			//Take inputs and store value for use
			int inputs;
			cin >> inputs;

			//Define variables for cell numbers, lengths, diffusion and source values
			int M;	
			int N;
			double Lx;
			double Ly;
			double alpha;
			double source;
			int filenumbers;

					//If user selected to use custom inputs, prompt and receive values for each
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

							//Number of files saved will correspond to the number of data points printed to vtk

							cout << endl << endl;
							cout << "Number of files to save: ";
							cin >> filenumbers;
						
							cout << endl << endl;
						break;

						default:
							//If any other value was input, use the default values hard-coded here
							M = 10;
							N = 10;
							Lx = 1;
							Ly = 2;
							alpha = 0.1;
							source = 0.1;
							filenumbers = 50;
						break;
					}

	//Create variables to store the resolutions in both directions, based on inputs

	double delx = Lx/M;
	double dely = Ly/N;


	vector<volume> volume_array;	//Create vector of volume class objects
	vector<face> face_array;		//Create vector of face class objects
	vector<double> Temp;			//Create vector of doubles to store Temperatures

	volume_array.reserve(M*N);				//Reserve space for CV's equal to x-cells times y-cells
	Temp.reserve(M*N);						//Reserve space for temperatures equal to number of CV's
	face_array.reserve(2*(M+1)*(N+1));		//Reserve space for enough faces correspondent to number of cells

	//Additional variables for face and cell numbering operations
	int k;
	int f;
	int p;


	//Loop through all cells with two loops
	for(int j = 0;j<N;j++){
		for(int i = 0;i<M;i++){
			k = ((i*N)+j);			//Setting cell value from 0 to M*N
			
			//Initializing values for control volumes
			volume_array[k].ID = k;
			volume_array[k].midpointx = (i*delx) + (delx/2);
			volume_array[k].midpointy = (j*dely) + (dely/2);

			//Setting normals for each face, which North and East positive
			volume_array[k].n_N = 1;
			volume_array[k].n_S = -1;
			volume_array[k].n_E = 1;
			volume_array[k].n_W = -1;


			//Setting array volume based on resolution
			volume_array[k].V = delx*dely;


			//Setting whether the control volume is on a boundary
			if(j == 0)			volume_array[k].bound_flag = 1;
			else if(j == N-1)	volume_array[k].bound_flag = 1;

			if(i== 0)			volume_array[k].bound_flag = 1;
			else if(i == M-1)	volume_array[k].bound_flag = 1;
		}
	}


	//Create variable for the boundary condition input and prompt for the input
	int boundarycondition;

	cout << endl << endl;
	cout << "Enter 1 for hard-coded boundary conditions, or anything else for the default periodic boundary conditions:  ";
	cin >> boundarycondition;


	//Creating and opening file for Face ID debugging
	ofstream Fout;
	Fout.open("faces.txt");





	switch(boundarycondition){
		//If boundary condition is set for hard-coding, following definitions are used 
		case 1:
			
			//Define east and west face IDs for non-periodic
			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);							//Face setting values equal with an offset for additional face per row
					face_array[k].ID = k;						
					face_array[k].A = dely;						//Setting face area to the resolution in the y-direction

					if(i==0) face_array[k].Neigh_1_ID = -1;		//If first face of row, set no neighbor to left
					else 	 face_array[k].Neigh_1_ID = k-1-j;

					if(i==M) face_array[k].Neigh_2_ID = -3;		//If last face of row, set no neighbor to right
					else 	 face_array[k].Neigh_2_ID = k-j;

				}
			}
			
			//Define east and west face IDs for non-periodic
			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);										//Setting same k value as before
					f = k+((M+1)*(N));										//F value to continue sequential ordering based on k value
					face_array[f].ID = f;
					face_array[f].A = delx;									//Setting face area to the resolution in the x-direction


					cout << k << " " << f << endl;							//Debugging line to see how k and f sequentially order

					if((j==0)&&(i!=M))	face_array[f].Neigh_1_ID = -2;		//If first face of column and not last member of the row, set no neighbor below
					else 		face_array[f].Neigh_1_ID =  k-M;			

					if((j==(N-1))&&(i!=0)) face_array[f].Neigh_2_ID = -4;	//If last face of column and not first member of the row, set no neighbor above
					else 	     face_array[f].Neigh_2_ID = k;
				}
			}


			//Print each face's ID and its neighbors to file for debugging and verification

			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);

					Fout << "Face ID: " << face_array[k].ID << "  N1: " << face_array[k].Neigh_1_ID << "  N2: " << face_array[k].Neigh_2_ID << endl;
				}
			}

			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);
					f = k+((M+1)*(N));

					Fout << "Face ID: " << face_array[f].ID << "  N1: " << face_array[f].Neigh_1_ID << "  N2: " << face_array[f].Neigh_2_ID << endl;

				}
			}

			Fout.close();

		break;




		
		//If any other input for BC's, use periodic boundary conditions outlines below
		default:
			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);							//Face setting values equal with an offset for additional face per row
					face_array[k].ID = k;						
					face_array[k].A = dely;						//Setting face area to the resolution in the y-direction

					if(i==0) face_array[k].Neigh_1_ID = (M-1)+(M*j);		//If first face of row, set end of row as neighbor
					else 	 face_array[k].Neigh_1_ID = k-1-j;

					if(i==M) face_array[k].Neigh_2_ID = (j*M);		//If last face of row, set first face of row as neighbor
					else 	 face_array[k].Neigh_2_ID = k-j;

				}
			}
			
			//Define north and south face IDs for non-periodic
			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);										//Setting same k value as before
					f = k+((M+1)*(N));										//F value to continue sequential ordering based on k value
					face_array[f].ID = f;
					face_array[f].A = delx;									//Setting face area to the resolution in the x-direction


					cout << k << " " << f << endl;							//Debugging line to see how k and f sequentially order

					if((j==0)&&(i!=M))	face_array[f].Neigh_1_ID = ((N-1)*N)+i;		//If first face of column and not last member of the row, last of column as neighbor
					else 		face_array[f].Neigh_1_ID =  k-M;			

					if((j==(N-1))&&(i!=0)) face_array[f].Neigh_2_ID = i-1;	//If last face of column and not first member of the row, first of column as neighbor
					else 	     face_array[f].Neigh_2_ID = k;
				}
			}


			//Print each face's ID and its neighbors to file for debugging and verification

			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);

					Fout << "Face ID: " << face_array[k].ID << "  N1: " << face_array[k].Neigh_1_ID << "  N2: " << face_array[k].Neigh_2_ID << endl;
				}
			}

			for(int j = 0;j<N;j++){
				for(int i = 0;i<=M;i++){
					k = ((j*(M+1))+i);
					f = k+((M+1)*(N));

					Fout << "Face ID: " << face_array[f].ID << "  N1: " << face_array[f].Neigh_1_ID << "  N2: " << face_array[f].Neigh_2_ID << endl;

				}
			}

			Fout.close();

		break;
	}






	

	//Establish timer variable for use in input of total runtime
	double timer;

	cout << endl << endl;
	cout << "Enter runtime for the simulation in seconds:  ";
	cin >> timer;

	//Set another variable for the timer of the length of the runtime
		clock_t timers;
		timers = clock();



	//Establishing timestep based on stabiity

	double delt = 0.1*(0.25*delx*dely);


	//Setting variables for initial condition
	double denom = sqrt(pow((Lx/2),2)+pow((Lx/2),2));
	double numa;
	double numb;

	//Create and open files for Temperature and Initial Conditions
	ofstream Tout;
	ofstream Iout;
  	Iout.open ("Data/Output.txt");
  	
  	//Calculate initial conditions as prescribed and print to file
	for(int j = 0;j<N;j++){
		for(int i = 0;i<M;i++){
			k = ((j*N)+i);

			numa = (volume_array[k].midpointx - (Lx/2));
			numb = (volume_array[k].midpointy - (Ly/2));

			numa = pow(numa,2);
			numb = pow(numb,2);

			Temp[k] = (numa + numb)/denom;

			Iout << Temp[k] << " ";
		}

		Iout << endl;
	}




	//Set counter variables for use in logging outputs to files
	int logcount = 0;
	int logger = timer/(filenumbers*delt);
	int storagecounter = 0;


	//Define variables for flux calculations
	double fluxN;
	double fluxS;
	double fluxE;
	double fluxW;
	double fluxsum;

	//Set strings for output vtk file names
	string s1 = "Data/data.";
	string s2 = ".vtk";
	string sint;
	string filename;


	//Loop for solution
	for(double t = delt;t<timer;t=t+delt){

			//If using hard-coded Boundary Conditions
				if(boundarycondition == 1){


					//Loop through flux calculations
					for(int j = 0;j<N;j++){
						for(int i = 0;i<=M;i++){
							k = ((j*(M+1))+i);
							f = k+((M+1)*(N));
 							
 							//If on right boundary (-1) use one-third source, on left (-3) use defined temperature
							if(face_array[k].Neigh_1_ID == -1)			face_array[k].flux = -abs(source/3);
							else if(face_array[k].Neigh_2_ID == -3) 	face_array[k].flux = ((dely*alpha)*(-1-(Temp[face_array[k].Neigh_1_ID])/2));
							//Otherwise use simple flux calculation
							else 										face_array[k].flux = face_array[k].cal_dif_flux_on_face(Temp,alpha);
 							
 							//If on bottom boundary (-2) use temperature, on top(-4) use defined flux of zero
							if(face_array[f].Neigh_1_ID == -2)			face_array[f].flux = ((delx*alpha)*((Temp[face_array[f].Neigh_2_ID]-1)/2));	
							else if(face_array[f].Neigh_2_ID == -4) 	face_array[f].flux = 0;
							//Otherwise use simple flux calculation
							else 										face_array[f].flux = face_array[f].cal_dif_flux_on_face(Temp,alpha);

						}
					}


					//Loop through cells and calculate sum of the flux and new temperature
					for(int j = 0;j<N;j++){
						for(int i = 0;i<M;i++){
							p = (j*N)+i;
							k = (j*(M))+i;
							f = k+((M+1)*(N));

						fluxN = face_array[f+M].flux;
						fluxS = face_array[f].flux;
						fluxE = face_array[k+1+j].flux;
						fluxW = face_array[k+j].flux;

						volume_array[p].fluxsum = volume_array[p].cal_flux_sum(fluxN,fluxS,fluxE,fluxW);
						Temp[p]= Temp[p] + ((delt/(delx*dely))*volume_array[p].fluxsum) + (delt*source);

						}
					}
				}


				//If periodic boundary conditions
				else{

					//Simple flux calculations all around
					for(int j = 0;j<N;j++){
						for(int i = 0;i<=M;i++){
							k = ((j*(M+1))+i);
							f = k+((M+1)*(N));
 							
 							face_array[k].flux = face_array[k].cal_dif_flux_on_face(Temp,alpha);
							face_array[f].flux = face_array[f].cal_dif_flux_on_face(Temp,alpha);

						}
					}

					//Loop through cells calculating flux sum and new temperature at each cell
					for(int j = 0;j<N;j++){
						for(int i = 0;i<M;i++){
							p = (j*N)+i;
							k = (j*(M))+i;
							f = k+((M+1)*(N));

						fluxN = face_array[f+M].flux;
						fluxS = face_array[f].flux;
						fluxE = face_array[k+1+j].flux;
						fluxW = face_array[k+j].flux;

						volume_array[p].fluxsum = volume_array[p].cal_flux_sum(fluxN,fluxS,fluxE,fluxW);
						Temp[p]= Temp[p] + ((delt/(delx*dely))*volume_array[p].fluxsum) + (delt*source);

						}
					}
				}


		//Check if the number of runs loope through equals the interval at which files should be saved
		if(logcount == logger){
			//If true, print the results to file 
			Iout << endl << endl;
			Iout << "Time: " << t << endl;
			Iout << "Data:" << endl;

			//Convert counter to string for use in file name
			sint = to_string(storagecounter);
			
			storagecounter = storagecounter + 1;

			filename = s1+sint+s2;

			//Create new file for each saving interval for vtk usage

			//File name is the concatenated strings
			Tout.open(filename.c_str());

			//Specific vtk file format requirements
			Tout << "# vtk DataFile Version 2.0" << endl;
			Tout << "Data for 2D FVM Solution" << endl;
			Tout << "ASCII" << endl; 
			Tout << "DATASET STRUCTURED_GRID" << endl;
			Tout << "DIMENSIONS " << M << " " << N  << " 1" << endl;
			Tout << "POINTS "<< (M*N) <<" float" << endl;

			//Print series of points as the locations for each cells with z coordinate of zero
			for(int j = N-1 ;j>=0;j--){
				for(int i = M-1;i>=0;i--){
					k = ((j*N)+i);
					Tout << volume_array[k].midpointx << " " << volume_array[k].midpointy << " 0"  << endl;
				}
			}

			//More formatting for vtk for temperature data
			Tout << "POINT_DATA " << (M*N) << endl;
			Tout << "SCALARS temperatures float" << endl;
			Tout << "LOOKUP_TABLE default" << endl;

			//Print temperature data to both files
			for(int j = N-1 ;j>=0;j--){
				for(int i = M-1;i>=0;i--){
					k = ((j*N)+i);
					Tout << Temp[k] << endl;
					Iout << Temp[k] << " ";
				}
				Iout << endl;
			}
			
			//Reset log counter
			logcount = 0;

			//Close current vtk file
			Tout.close();

			//Print to terminal current progress
			cout << "Progress: (" << storagecounter << "/" << filenumbers << ")" << endl;
		}
		else logcount = logcount+1;

	}

	//Finish printing of results when loop is finished
	Iout << endl;
	Iout << "Final Results: " << endl;

	for(int j = N-1 ;j>=0;j--){
		for(int i = M-1;i>=0;i--){
			k = ((j*N)+i);

			Iout << Temp[k] << " ";
		}

		Iout << endl;
	}

	Iout.close();

	//Calculate total time taken and print to terminal
	timers = clock() - timers;

	cout << endl << endl;
			
	cout << "\nTime elapsed: "<< (float)timers/CLOCKS_PER_SEC << " seconds" << endl << endl;

	cout << endl << endl;

	return 0;
}