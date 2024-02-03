/******************************************************************************

C subroutine for a Navier Slip Law
Each wall is defined as a structure 
The wall characteristics are stored in an array of structures afterwards

*******************************************************************************/

#include "udf.h"

//Maximum number of walls
#define NUMMAXWALLS 10

/******************************************************************************

Wall structure

*******************************************************************************/

// Structure of a wall
struct WallInfo {
    int ID;                   // Zone ID of the wall
    int MovType;              // Movement type
    double Vector_1[3];       // Vectors of the wall movement
    double Vector_2[3];
    double BetaAdim;          // Non-dimensional slip coefficient
    double Relax_factor;      // Under-relaxing factor
    double href;              // Length of reference
};

// Global array to store wall information
struct WallInfo Walls[NUMMAXWALLS];


/******************************************************************************

Void function to fill the array of structures

*******************************************************************************/

void fill_array_structure()
{
    FILE *file;
    char filename[] = "C:/Users/flemarchand/Documents/Prueba_UDF_slip_navier/udf_informe/udf_moving_wall_slip/txt_file_walls.txt";
    char line[500];

    // Open the file
    file = fopen(filename, "r");
    if (file == NULL) {
        //Message("Error: File can't be opened.\n");
        return;  // Exit if file cannot be opened
    }

    // Read the file and assign the values to the array of structures
	int i = 0;
	while (fgets(line, sizeof(line), file) != NULL)  // fgets to read the whole line
	{
		// Skip empty lines or lines that start with a comment character (header)
		if (line[0] == '\n') continue;
		
        // Parse the line as a wall data entry
		if(sscanf(line, "(%d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf)",
				   &Walls[i].ID, &Walls[i].MovType,
				   &Walls[i].Vector_1[0], &Walls[i].Vector_1[1], &Walls[i].Vector_1[2],
				   &Walls[i].Vector_2[0], &Walls[i].Vector_2[1], &Walls[i].Vector_2[2],
				   &Walls[i].BetaAdim, &Walls[i].Relax_factor, &Walls[i].href) == 11){

			// Check if the movement is set correctly
			switch (Walls[i].MovType) {
				case 0:
					//Message("Translational movement for Wall ID: %d\n", Walls[i].ID);
					break;
				case 1:
					//Message("Rotational movement for Wall ID: %d\n", Walls[i].ID);
					break;
				default:
					//Message("Error: Movement type not set correctly for Wall ID: %d\n", Walls[i].ID);
					break;
			}

			i++; // Increment the index only if the line was successfully parsed
			if (i >= NUMMAXWALLS) break; // Prevent reading beyond the array bounds
		}else{
			// Parsing failed
			//Message("Error parsing line: %s\n", line);
		}
	}
	fclose(file);
}


/****************************************************************************************************************************************

Store in F_UDMI the calculated velocities and initialise the betas of the walls using the on_demand macro, 4 F_UDMI in total at this step

****************************************************************************************************************************************/

DEFINE_ON_DEMAND(on_demand_calc)
{
   #if RP_NODE
   //Use the void function to access the structure information
   fill_array_structure();
   
   Domain *domain; /* declare domain pointer since it is not passed as an argument to the DEFINE macro */
   
   Thread *t;
   face_t f;
   domain = Get_Domain(1);  /* Get the domain using Ansys Fluent utility */
   
   //Cell thread and viscosity definition needed to access the viscosity values in contour
   cell_t c0;
   Thread *t0;
   real mu_c;
   
   //Local variables needed to calculate the velocities if needed (case of a rotational velocity)
   real vw[3]={0.0,0.0,0.0};
   real omega[3]={0.0,0.0,0.0};
   real x0[3]={0.0,0.0,0.0};
   real x[ND_ND];
   real a[ND_ND];
   
   for(int i=0;i<NUMMAXWALLS;i++){
	   //Check that only walls ID are used
	   if(Walls[i].ID!=0){
			t=Lookup_Thread(domain,Walls[i].ID);
			begin_f_loop(f,t)
			{
				//Velocity calculation and storage in F_UDMI
				//Case of a translational velocity, assign to the F_UDMI the values of the Vector 1
				if(Walls[i].MovType==0){
					F_UDMI(f,t,0)=Walls[i].Vector_1[0];
					F_UDMI(f,t,1)=Walls[i].Vector_1[1];
					F_UDMI(f,t,2)=Walls[i].Vector_1[2];
	
				}else{
				
					//Case of a rotational velocity, calculation of the velocity components in function of the rotational velocity
					for(int j=0;j<3;j++){
						omega[j]=Walls[i].Vector_1[j];
						x0[j]=Walls[i].Vector_2[j];
					}
					F_CENTROID(x, f, t);
					NV_VV(a, =, x, -, x0);
					NV_CROSS(vw, omega, a);
					F_UDMI(f,t,0)=vw[0];
					F_UDMI(f,t,1)=vw[1];
					F_UDMI(f,t,2)=vw[2];
		
				}
				//Beta calculation and storage in F_UDMI
				//Neighbour cell thread
				c0=F_C0(f,t);
				t0=THREAD_T0(t);
				mu_c=C_MU_L(c0,t0);
				F_UDMI(f,t,3)=((Walls[i].BetaAdim)*mu_c)/Walls[i].href;
			}
			end_f_loop(f,t)
	   }
	}
	#endif
} 

/************************************************************************************************************************************

Define_PROFILE that calculates the shear stress given a slip coefficient and a slip velocity for either a moving or a stationary wall

************************************************************************************************************************************/

#define UDM_POS_SHEARSTRESS 4
#define F_ShearStress(f_foo,th_foo,d_foo) F_UDMI(f_foo,th_foo,UDM_POS_SHEARSTRESS+d_foo)


//Function declaration
void GetShearStress(Thread*,int,int);

//Shear stress calculation
void GetShearStress(Thread* thread,int direction,int position)
{
	#if RP_NODE
	//Use the void function to access the structure information
    fill_array_structure();
	
	//Variables that will be employed to calculate the shear stress
	//Vectors used to calculate the rotational velocity
	real x[ND_ND];
	real x0[3]={0.0,0.0,0.0};
	real a[ND_ND];
	real omega[3]={0.0,0.0,0.0};
	
	//Vectors used to calculate the slip velocity
	real NV_VEC(vf);
	real NV_VEC(vw);
	real NV_VEC(vslip);
	real NV_VEC(vrel);
	real NV_VEC(vreln);
	real NV_VEC(Area);
	real NV_VEC(An);
	
	//Vectors used to calculate the shear stress
	real NV_VEC(tau);
	real NV_VEC(tau2);
	real beta;
	
	//Definition of a zero vector to use the NV_V_VS macro 
	real nulle[3]={0.0,0.0,0.0};
	
	//Cell thread needed to access the viscosity values in contour, definition of the viscosity variable
	face_t f;
	cell_t c0;
	Thread *t0;
	real mu_c;
	
	//Access the thread ID and use the corresponding values stored in the structure
	for(int i=0;i<NUMMAXWALLS;i++){
		if((THREAD_ID(thread)==Walls[i].ID)&&(Walls[i].ID!=0)){
			//f_loop and initialisation using the F_UDMI
			begin_f_loop(f, thread){
				vw[0]=F_UDMI(f,thread,0);
				vw[1]=F_UDMI(f,thread,1);
				vw[2]=F_UDMI(f,thread,2);
				beta=F_UDMI(f,thread,3);
			}
			end_f_loop(f,thread)
			
			
			//f_loop and slip velocity calculation
			begin_f_loop(f, thread){
				
				//Area vector to get the normal vector components
				F_AREA(Area,f,thread);
				
				//Normalisation
				NV_V_VS(An,=,nulle,+,Area,*,(1/NV_MAG(Area)));
				
				/* Get the velocities of the fluid in the boundary */
				vf[0]=F_U(f,thread);
				vf[1]=F_V(f,thread);
				vf[2]=F_W(f,thread);
				
				//Calculation of the tangential slip velocity in function of the movement type
				if(Walls[i].MovType==0){
					vw[0]=F_UDMI(f,thread,0);
					vw[1]=F_UDMI(f,thread,1);
					vw[2]=F_UDMI(f,thread,2);
					NV_VV(vrel, =, vf, -, vw);
					NV_V_VS(vreln, =,nulle,+,An,*,(NV_DOT(vrel, An)));
					NV_VV(vslip, =, vrel, -, vreln);
				}else if(Walls[i].MovType==1){
					//Case of a set rotational velocity
					for(int j=0;j<3;j++){
						omega[j]=Walls[i].Vector_1[j];
						x0[j]=Walls[i].Vector_2[j];
					}
					F_CENTROID(x, f, thread);
					NV_VV(a, =, x, -, x0);
					NV_CROSS(vw, omega, a);
					NV_VV(vrel, =, vf, -, vw);
					NV_V_VS(vreln, =,nulle,+,An,*,(NV_DOT(vrel, An)));
					NV_VV(vslip, =, vrel, -, vreln);
				}
				
				//Shear stress calculation
				t0=THREAD_T0(thread);
				c0=F_C0(f,thread);
				mu_c=C_MU_L(c0,t0);
				beta=((Walls[i].BetaAdim)*mu_c)/Walls[i].href;
				NV_V_VS(tau, =,nulle,+,vslip,*,beta);
				
				//under-relaxation
				tau2[direction]=(1.0-Walls[i].Relax_factor)*F_ShearStress(f,thread,direction)+Walls[i].Relax_factor*tau[direction];
				F_ShearStress(f,thread,direction)=tau2[direction];
				F_PROFILE(f, thread, position) = tau2[direction];
			}
			end_f_loop(f,thread)
			
		}
	}
	//End of the profile definition
	#endif
}

//Define profile applied on each component
DEFINE_PROFILE(wall_slip_navier_x, thread, position)
{
	GetShearStress(thread,0,position);
}

DEFINE_PROFILE(wall_slip_navier_y, thread, position)
{
	GetShearStress(thread,1,position);
}

DEFINE_PROFILE(wall_slip_navier_z, thread, position)
{
	GetShearStress(thread,2,position);
}

/******************************************************************************

End of the subroutine

*******************************************************************************/
