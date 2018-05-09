/* Source File: Computational Physics 2018
 * liesegang_growth.cpp
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#include "liesegang.h"


/*'dimensionalized denotes
 * whether 'normal parameters' are used
 * or the advanced ones'
 *
 * This method implements the Coupled Set of differential equations
 * that describe the chemical process described by the
 * Liesegang phenomenon*/
int Liesegang::LiesegangGrowth(int mode)
{

	/*Mode for Equations WITH dimensions*/
	if(mode == 0)
	{

		double inverse_spatial = (1/double(pow(delta_X_,2)));

		/*Calculation of alpha factors*/
		alpha_a_ = D_a_ * delta_T_ * inverse_spatial;
		alpha_b_ = D_b_ * delta_T_ * inverse_spatial;
		alpha_c_ = D_c_ * delta_T_ * inverse_spatial;

		print("CFL: " , alpha_a_, "\n");

		/*Time Step Loop*/
		for(double t = 1*(delta_T_); t <= double(DEF_EXC_TIME); t = t + delta_T_)
		{
			/*Initialies every new time step with the template vector*/
			tube_.push_back(tube_template_);


			/*Restricts Output to later times*/
			if(t >= (double(DEF_EXC_TIME)-100*delta_T_))
			{
				a_concentration_hist << tube_.back().at(0).a_concentration << ",";
				b_concentration_hist << tube_.back().at(0).b_concentration << ",";
				c_concentration_hist <<tube_.back().at(0).c_concentration << ",";
				s_concentration_hist <<tube_.back().at(0).s_concentration << ",";
			}

			/*Spatial Loop - Iterates over Tube Units
			 * -> CONSIDERS BOUNDARY CONDITIONS
			 * -> TODO: VERIFY BOUNDARY CONDITIONS REALIZED AS LOOP BOUNDARIES*/
			for(int idx = 1; idx < DEF_GRID_N-1; idx++)
			{






				/*CONCENTRATION A */
				tube_.back().at(idx).a_concentration =
						tube_.at(0).at(idx).a_concentration
						+ alpha_a_* (tube_.at(0).at(idx+1).a_concentration
								- 2 * tube_.at(0).at(idx).a_concentration
								+ tube_.at(0).at(idx-1).a_concentration)
								-R_*delta_T_*tube_.at(0).at(idx).a_concentration*tube_.at(0).at(idx).b_concentration	;

				/*CONCENTRATION B */
				tube_.back().at(idx).b_concentration =
						tube_.at(0).at(idx).b_concentration
						+ alpha_b_* (tube_.at(0).at(idx+1).b_concentration
								- 2 * tube_.at(0).at(idx).b_concentration
								+ tube_.at(0).at(idx-1).b_concentration)
								-R_*delta_T_*tube_.at(0).at(idx).a_concentration*tube_.at(0).at(idx).b_concentration;







				tube_.back().at(idx).c_concentration =
						tube_.at(0).at(idx).c_concentration
						+ alpha_c_* (tube_.at(0).at(idx+1).c_concentration
								- 2 * tube_.at(0).at(idx).c_concentration
								+ tube_.at(0).at(idx-1).c_concentration)
								+ R_*delta_T_*tube_.at(0).at(idx).a_concentration*tube_.at(0).at(idx).b_concentration
								- N_one_*delta_T_* Heaviside(tube_.at(0).at(idx).c_concentration)*
								(tube_.at(0).at(idx).c_concentration*tube_.at(0).at(idx).c_concentration)
								- N_two_*delta_T_*tube_.at(0).at(idx).c_concentration * tube_.at(0).at(idx).s_concentration;



				/*CONCENTRATION S -> Boundary Conditions not Important*/
				tube_.back().at(idx).s_concentration =
						tube_.at(0).at(idx).s_concentration
						+ N_one_ *Heaviside(tube_.at(0).at(idx).c_concentration)*delta_T_*
						(tube_.at(0).at(idx).c_concentration*tube_.at(0).at(idx).c_concentration)
						+ N_two_ * tube_.at(0).at(idx).c_concentration* tube_.at(0).at(idx).s_concentration*delta_T_;


				/*Restricts Output to later times*/


				/*Creates File Ouput for Data Analysis*/
				if(t >= (double(DEF_EXC_TIME)-100*delta_T_))
				{
					a_concentration_hist << tube_.back().at(idx).a_concentration << ",";
					b_concentration_hist << tube_.back().at(idx).b_concentration << ",";
					c_concentration_hist <<tube_.back().at(idx).c_concentration << ",";
					s_concentration_hist <<tube_.back().at(idx).s_concentration << ",";
				}

				tube_.back().at(idx).time_step = I_time_step_;



			}

			/*Restricts Output to later times*/

			if(t >= (double(DEF_EXC_TIME)-100*delta_T_))
			{
				a_concentration_hist << tube_.back().back().a_concentration ;
				b_concentration_hist << tube_.back().back().b_concentration ;
				c_concentration_hist <<tube_.back().back().c_concentration;
				s_concentration_hist <<tube_.back().back().s_concentration ;

				a_concentration_hist << "\n";
				b_concentration_hist << "\n";
				c_concentration_hist << "\n";
				s_concentration_hist << "\n";

			}

			/*Time Step Increment*/
			I_time_step_ = I_time_step_+delta_T_;


			/*Erases not needed vector element @I_time_step-2 -> first element
			 * This saves a looooooooot of memory -> foobar*/
			tube_.erase(tube_.begin());

		}

	}


	/*Implements use of diffusion functional*/
	else if (mode == 1)
	{




		double inverse_spatial = (1/double(pow(delta_X_,2)));


		/*RE-Calculation of Diffutions Functionals (former constants)
		 * -> USE of s-values from former time step
		 * -> alpha factors are now used for another 'concentration'
		 * -> alpha's do not include delta_t anymore*/


		s_zero_ = double(DEF_S_ZERO);





		print("CFL: " , alpha_a_, "\n");

		/*Time Step Loop*/
		for(double t = 1*(delta_T_); t <= double(DEF_EXC_TIME); t = t + delta_T_)
		{
			/*Initialies every new time step with the template vector*/
			tube_.push_back(tube_template_);


			/*Restricts Output to later times*/


			if(t >= (double(DEF_EXC_TIME)*299/300))
			{
				a_concentration_hist << tube_.back().at(0).a_concentration << ",";
				b_concentration_hist << tube_.back().at(0).b_concentration << ",";
				c_concentration_hist <<tube_.back().at(0).c_concentration << ",";
				s_concentration_hist <<tube_.back().at(0).s_concentration << ",";

			}
			/*Spatial Loop - Iterates over Tube Units
			 * mode 1 ->
			 * changes equation system*/
			for(int idx = 1; idx < DEF_GRID_N-1; idx++)
			{



				/*RE-Calculation of alpha-factors
				 * -> used for P_block now
				 * -> do NOT include timestep anymoe in this step*/

				alpha_a_ = (1/1+(tube_.at(0).at(idx).s_concentration/s_zero_)) * D_a_
						* inverse_spatial;
				alpha_b_ = (1/1+(tube_.at(0).at(idx).s_concentration/s_zero_)) * D_b_
						* inverse_spatial;
				alpha_c_ = (1/1+(tube_.at(0).at(idx).s_concentration/s_zero_)) * D_c_
						* inverse_spatial;

				/*------------------------------------------*/


				/*'CONCENTRATION' of respective P_blocks */
				/*A*/	tube_.back().at(idx).P_block_a =
						alpha_a_ * (tube_.at(0).at(idx).a_concentration - tube_.at(0).at(idx-1).a_concentration);

				/*B*/	tube_.back().at(idx).P_block_b=
						alpha_b_ * (tube_.at(0).at(idx).b_concentration - tube_.at(0).at(idx-1).b_concentration);

				/*C*/	tube_.back().at(idx).P_block_c =
						alpha_c_ * (tube_.at(0).at(idx).c_concentration - tube_.at(0).at(idx-1).c_concentration);


				/*TODO: Introduce Boundary Conditions to P_Blocks MAYBE*/






				/*CONCENTRATION A */
				tube_.back().at(idx).a_concentration =
						tube_.at(0).at(idx).a_concentration
						+ delta_T_ * inverse_spatial *
						(tube_.back().at(idx).P_block_a - tube_.back().at(idx-1).P_block_a)
						-R_*delta_T_*tube_.at(0).at(idx).a_concentration*tube_.at(0).at(idx).b_concentration	;

				/*CONCENTRATION B */
				tube_.back().at(idx).b_concentration =
						tube_.at(0).at(idx).b_concentration
						+ delta_T_ * inverse_spatial *
						(tube_.back().at(idx).P_block_b - tube_.back().at(idx-1).P_block_b)
						-R_*delta_T_*tube_.at(0).at(idx).a_concentration*tube_.at(0).at(idx).b_concentration;




				tube_.back().at(idx).c_concentration =
						tube_.at(0).at(idx).c_concentration
						+ delta_T_ * inverse_spatial *
						(tube_.back().at(idx).P_block_c - tube_.back().at(idx-1).P_block_c)
						+ R_*delta_T_*tube_.at(0).at(idx).a_concentration*tube_.at(0).at(idx).b_concentration
						- N_one_*delta_T_* Heaviside(tube_.at(0).at(idx).c_concentration)*
						(tube_.at(0).at(idx).c_concentration*tube_.at(0).at(idx).c_concentration)
						- N_two_*delta_T_*tube_.at(0).at(idx).c_concentration * tube_.at(0).at(idx).s_concentration;



				/*CONCENTRATION S -> Boundary Conditions not Important*/
				tube_.back().at(idx).s_concentration =
						tube_.at(0).at(idx).s_concentration
						+ N_one_ *Heaviside(tube_.at(0).at(idx).c_concentration)*delta_T_*
						(tube_.at(0).at(idx).c_concentration*tube_.at(0).at(idx).c_concentration)
						+ N_two_ * tube_.at(0).at(idx).c_concentration* tube_.at(0).at(idx).s_concentration*delta_T_;


				/*Restricts Output to later times*/


				/*Creates File Ouput for Data Analysis*/
				if(t >= (double(DEF_EXC_TIME)*3/4))
				{
					a_concentration_hist << tube_.back().at(idx).a_concentration << ",";
					b_concentration_hist << tube_.back().at(idx).b_concentration << ",";
					c_concentration_hist <<tube_.back().at(idx).c_concentration << ",";
					s_concentration_hist <<tube_.back().at(idx).s_concentration << ",";
				}

				tube_.back().at(idx).time_step = I_time_step_;


			}

			/*Restricts Output to later times*/

			if(t >= (double(DEF_EXC_TIME)*3/4))
			{
				a_concentration_hist << tube_.back().back().a_concentration ;
				b_concentration_hist << tube_.back().back().b_concentration ;
				c_concentration_hist <<tube_.back().back().c_concentration;
				s_concentration_hist <<tube_.back().back().s_concentration ;

				a_concentration_hist << "\n";
				b_concentration_hist << "\n";
				c_concentration_hist << "\n";
				s_concentration_hist << "\n";

			}

			/*Time Step Increment*/
			I_time_step_ = I_time_step_+delta_T_;


			/*Erases not needed vector element @I_time_step-2 -> first element
			 * This saves a looooooooot of memory -> foobar*/
			tube_.erase(tube_.begin());

		}

	}

	return 0;
}
