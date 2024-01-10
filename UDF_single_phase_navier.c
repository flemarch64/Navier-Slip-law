//Single phase version

#include "udf.h"
#include "sg.h"


DEFINE_PROFILE(wall_slip_navier_stationnary, thread, position)
{
	//UDF to calculate the shear stress by implementing the Navier Slip law at a stationary wall
	
	#if RP_NODE
	
	//Schemes variables to be employed
	//Under relaxation factor
	real alfa_relax=RP_Get_Real("alfa_relax");
	
	//Non-dimensional Slip coefficient for the liquid to be tested
	real beta1_adim=RP_Get_Real("beta1_adim");
	
	//beta dimensionalisation, with viscosity and href
	real mu=RP_Get_Real("mu");
	real href=RP_Get_Real("h_ref");
	//End of the schemes variables definition
	
	//Profile definition
	face_t f;
	real x[ND_ND];
	real NV_VEC(vf);
	real NV_VEC(vfn);
	real NV_VEC(Area);
	real NV_VEC(tau);
	real NV_VEC(tau2);
	real NV_VEC(vslip);
	real NV_VEC(vrel);
	real NV_VEC(vreln);
	real nulle[3]={0.0,0.0,0.0};
	real NV_VEC(An);
	real beta1;
	beta1=(beta1_adim*mu)/href;

	begin_f_loop(f, thread){
		//Area vector to get the normal vector components
		F_AREA(Area,f,thread);
		
		//Normalisation
		NV_V_VS(An,=,nulle,+,Area,*,(1/NV_MAG(Area)));
		
		/* Get the velocities of the fluid in the boundary */
		vf[0]=F_U(f,thread);
		vf[1]=F_V(f,thread);
		vf[2]=F_W(f,thread);
		
		//Fluid velocity
		NV_V_VS(vfn, =,nulle,+,An,*,(NV_DOT(vf, An)));
		NV_VV(vslip, =, vf, -, vfn);
;
		
		//Shear stress calculation
		NV_V_VS(tau, =,nulle,+,vslip,*,beta1);
		
		//Switch case to assign the shear stress components in the GUI, underelaxation
		switch(position){
			case 0:
				tau2[0]=(1.0-alfa_relax)*F_UDMI(f,thread,0)+alfa_relax*tau[0];
				F_UDMI(f,thread,0)=tau2[0];
				F_PROFILE(f, thread, position) = tau2[0];
				break;
			case 1:
				tau2[1]=(1.0-alfa_relax)*F_UDMI(f,thread,1)+alfa_relax*tau[1];
				F_UDMI(f,thread,1)=tau2[1];
				F_PROFILE(f, thread, position) = tau2[1];
				break;
			case 2:
				tau2[2]=(1.0-alfa_relax)*F_UDMI(f,thread,2)+alfa_relax*tau[2];
				F_UDMI(f,thread,2)=tau2[2];
				F_PROFILE(f, thread, position) = tau2[2];
				break;
		}
	}
	end_f_loop(f,thread)
	//End of the profile definition
	#endif
}


DEFINE_PROFILE(wall_slip_navier_moving, thread, position)
{
	#if RP_NODE
	
	//UDF to calculate the shear stress by implementing the Navier slip law for a moving wall
	
	//Schemes variables to be employed
	//Under relaxation factor
	real alfa_relax=RP_Get_Real("alfa_relax");
	
	//Non-dimensional Slip coefficient for the liquid to be tested
	real beta1_adim=RP_Get_Real("beta1_adim");
	
	//Definition of a wall movement, general case
	int movement_type;
	movement_type=RP_Get_Integer("movement_type");
	
	//Case of a translational velocity
	real translational_velocity=RP_Get_Real("translation_velocity_modulus");
	real vx_dir=RP_Get_Real("vec_trans_direction_x");
	real vy_dir=RP_Get_Real("vec_trans_direction_y");
	real vz_dir=RP_Get_Real("vec_trans_direction_z");
	real vector_direction[3]={vx_dir,vy_dir,vz_dir};
	
	//Case of a rotational velocity
	real rotational_velocity=RP_Get_Real("rotation_velocity_modulus");
	real vec_rot_x_axis=RP_Get_Real("vec_rot_x");
	real vec_rot_y_axis=RP_Get_Real("vec_rot_y");
	real vec_rot_z_axis=RP_Get_Real("vec_rot_z");
	real rot_axis[3]={vec_rot_x_axis,vec_rot_y_axis,vec_rot_z_axis};
	
	//beta dimensionalisation, with viscosity and href
	real mu=RP_Get_Real("mu");
	real href=RP_Get_Real("h_ref");
	//End of the schemes variables definition
	
	//Profile definition
	face_t f;
	real x[ND_ND];
	real omega[3]={0.0,0.0,0.0};
	real NV_VEC(vf);
	real NV_VEC(vfn);
	real vw[3]={0.0,0.0,0.0};
	real NV_VEC(Area);
	real NV_VEC(tau);
	real NV_VEC(tau2);
	real NV_VEC(vslip);
	real NV_VEC(vrel);
	real NV_VEC(vreln);
	real nulle[3]={0.0,0.0,0.0};
	real NV_VEC(An);
	real beta1;
	beta1=(beta1_adim*mu)/href;

	//If to take into account the moving wall in the zone ID
	begin_f_loop(f, thread){
	
		//Area vector to get the normal vector components
		F_AREA(Area,f,thread);
		
		//Normalisation
		NV_V_VS(An,=,nulle,+,Area,*,(1/NV_MAG(Area)));
		
		/* Get the velocities of the fluid in the boundary */
		vf[0]=F_U(f,thread);
		vf[1]=F_V(f,thread);
		vf[2]=F_W(f,thread);
		
		//Centroids calculation to calculate the wall velocity in case of a rotating wall
		F_CENTROID(x, f, thread);
		
		//Movement type, wall velocity calculation
		switch (movement_type){
			case 0:
				NV_V_VS(vw, =,nulle,+,vector_direction,*,translational_velocity);
				break;
			case 1:
				NV_V_VS(omega, =,nulle,+,rot_axis,*,rotational_velocity);
				NV_CROSS(vw, x, omega);
				break;
		}
		
		//Relative velocity (normal component)
		NV_VV(vrel, =, vf, -, vw);
		NV_V_VS(vreln, =,nulle,+,An,*,(NV_DOT(vrel, An)));
		NV_VV(vslip, =, vrel, -, vreln);

		//Shear stress calculation
		NV_V_VS(tau, =,nulle,+,vslip,*,beta1);
	
		//Switch case to assign the shear stress components in the GUI, underelaxation
		switch(position){
			case 0:
				tau2[0]=(1.0-alfa_relax)*F_UDMI(f,thread,0)+alfa_relax*tau[0];
				F_UDMI(f,thread,0)=tau2[0];
				F_PROFILE(f, thread, position) = tau2[0];
				break;
			case 1:
				tau2[1]=(1.0-alfa_relax)*F_UDMI(f,thread,1)+alfa_relax*tau[1];
				F_UDMI(f,thread,1)=tau2[1];
				F_PROFILE(f, thread, position) = tau2[1];
				break;
			case 2:
				tau2[2]=(1.0-alfa_relax)*F_UDMI(f,thread,2)+alfa_relax*tau[2];
				F_UDMI(f,thread,2)=tau2[2];
				F_PROFILE(f, thread, position) = tau2[2];
				break;
		}
	}
	end_f_loop(f,thread)
	#endif
}
