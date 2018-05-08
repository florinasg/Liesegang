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
int Liesegang::LiesegangGrowth(bool dimensionalized)
{

	/*Calculation of alpha factors*/
	alpha_a_ = D_a_ * delta_T_ * (1/double(pow(delta_X_,2)));
	alpha_b_ = D_b_ * delta_T_ * (1/double(pow(delta_X_,2)));
	alpha_c_ = D_c_ * delta_T_ * (1/double(pow(delta_X_,2)));


	/*Time Step Loop*/
	for(double t = 1*(delta_T_); t <= double(DEF_EXC_TIME); t = t * delta_T_)
	{
		/*Initialies every new time step with the template vector*/
		tube_.push_back(tube_template_);

		/*Spatial Loop - Iterates over Tube Units
		 * -> CONSIDERS BOUNDARY CONDITIONS*/
		for(int idx = 1; idx < DEF_GRID_N-1; idx++)
		{


			/*Special Case*/
			if(idx == 1)
			{
				/*CONCENTRATION A */
				tube_.back().at(idx).a_concentration =
						tube_.at(I_time_step_-1).at(idx).a_concentration
						+ alpha_a_* (tube_.at(I_time_step_-1).at(idx+1).a_concentration
								- 2 * tube_.at(I_time_step_-1).at(idx).a_concentration)
								-R_*tube_.at(I_time_step_-1).at(idx).a_concentration*tube_.at(I_time_step_-1).at(idx).b_concentration;

				/*CONCENTRATION B*/
				tube_.back().at(idx).b_concentration =
						tube_.at(I_time_step_-1).at(idx).b_concentration
						+ alpha_b_* (tube_.at(I_time_step_-1).at(idx+1).b_concentration
								- 2 * tube_.at(I_time_step_-1).at(idx).b_concentration)
								-R_*tube_.at(I_time_step_-1).at(idx).a_concentration*tube_.at(I_time_step_-1).at(idx).b_concentration;;

				/*-------------------------------------------------------------------------*/

				/*CONCENTRATION C*/
				tube_.back().at(idx).c_concentration =
						tube_.at(I_time_step_-1).at(idx).c_concentration
						+ alpha_c_* (tube_.at(I_time_step_-1).at(idx+1).c_concentration
								- 2 * tube_.at(I_time_step_-1).at(idx).c_concentration)
								- N_one_* Heaviside(tube_.at(I_time_step_-1).at(idx).c_concentration)*
								(double(pow(tube_.at(I_time_step_-1).at(idx).c_concentration,2)))
								- N_two_*tube_.at(I_time_step_-1).at(idx).c_concentration * tube_.at(I_time_step_-1).at(idx).s_concentration;

			}


			if(idx == DEF_GRID_N-1)
			{
				/*CONCENTRATION A */
				tube_.back().at(idx).a_concentration =
						tube_.at(I_time_step_-1).at(idx).a_concentration
						+ alpha_a_*
						(- 2 * tube_.at(I_time_step_-1).at(idx).a_concentration
								+ tube_.at(I_time_step_-1).at(idx-1).a_concentration)
								-R_*tube_.at(I_time_step_-1).at(idx).a_concentration*tube_.at(I_time_step_-1).at(idx).b_concentration;

				/*CONCENTRATION B */
				tube_.back().at(idx).b_concentration =
						tube_.at(I_time_step_-1).at(idx).b_concentration
						+ alpha_b_*
						(- 2 * tube_.at(I_time_step_-1).at(idx).b_concentration
								+ tube_.at(I_time_step_-1).at(idx-1).b_concentration)
								-R_*tube_.at(I_time_step_-1).at(idx).a_concentration*tube_.at(I_time_step_-1).at(idx).b_concentration;;

				/*-------------------------------------------------------------------------*/

				/*CONCENTRATION C*/
				tube_.back().at(idx).c_concentration =
						tube_.at(I_time_step_-1).at(idx).c_concentration
						+ alpha_c_*
						(- 2 * tube_.at(I_time_step_-1).at(idx).c_concentration
								+ tube_.at(I_time_step_-1).at(idx-1).c_concentration)
								- N_one_* Heaviside(tube_.at(I_time_step_-1).at(idx).c_concentration)*
								(double(pow(tube_.at(I_time_step_-1).at(idx).c_concentration,2)))
								- N_two_*tube_.at(I_time_step_-1).at(idx).c_concentration * tube_.at(I_time_step_-1).at(idx).s_concentration;



			}

			else
			{

				/*CONCENTRATION A */
				tube_.back().at(idx).a_concentration =
						tube_.at(I_time_step_-1).at(idx).a_concentration
						+ alpha_a_* (tube_.at(I_time_step_-1).at(idx+1).a_concentration
								- 2 * tube_.at(I_time_step_-1).at(idx).a_concentration
								+ tube_.at(I_time_step_-1).at(idx-1).a_concentration)
								-R_*tube_.at(I_time_step_-1).at(idx).a_concentration*tube_.at(I_time_step_-1).at(idx).b_concentration	;

				/*CONCENTRATION B */
				tube_.back().at(idx).a_concentration =
						tube_.at(I_time_step_-1).at(idx).b_concentration
						+ alpha_b_* (tube_.at(I_time_step_-1).at(idx+1).b_concentration
								- 2 * tube_.at(I_time_step_-1).at(idx).b_concentration
								+ tube_.at(I_time_step_-1).at(idx-1).b_concentration)
								-R_*tube_.at(I_time_step_-1).at(idx).a_concentration*tube_.at(I_time_step_-1).at(idx).b_concentration;;

				/*-------------------------------------------------------------------------*/

				/*CONCENTRATION C*/
				tube_.back().at(idx).c_concentration =
						tube_.at(I_time_step_-1).at(idx).c_concentration
						+ alpha_c_* (tube_.at(I_time_step_-1).at(idx+1).c_concentration
								- 2 * tube_.at(I_time_step_-1).at(idx).c_concentration
								+ tube_.at(I_time_step_-1).at(idx-1).c_concentration)
								- N_one_* Heaviside(tube_.at(I_time_step_-1).at(idx).c_concentration)*
								(double(pow(tube_.at(I_time_step_-1).at(idx).c_concentration,2)))
								- N_two_*tube_.at(I_time_step_-1).at(idx).c_concentration * tube_.at(I_time_step_-1).at(idx).s_concentration;






			}








		}

		/*Time Step Increment*/
		I_time_step_ = I_time_step_+1;

	}


	return 0;
}
