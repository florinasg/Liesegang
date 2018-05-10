/* Source File: Computational Physics 2018
 * IinitializeTube.cpp
 *
 *  Created on: 08.05.2018
 *      Author: Sara Gasparini
 */


#include "liesegang.h"
#include <iostream>

int Liesegang::InitializeTube()
{

	/*Calculate DELTA_X as distance between grid-points*/
	delta_X_ = (1.0 /
			double(DEF_GRID_N-1));


	delta_T_ = double(DEF_DELTA_T);


	I_time_step_ = 1; /*Initialization */

	/*Container for Initial System ->
	 * serves also as template for every subsequent system vector*/
	tube_.push_back(std::vector<TubeSpaceUnit_>());



	/*Initializes Tube for Time Step 0*/
	for(int i  = 0; i < DEF_GRID_N; i = i + 1)
	{
		tube_.back().push_back(TubeSpaceUnit_());
		tube_.back().back().position = i*delta_X_; /*Stores position of tube unit*/
		tube_.back().back().time_step = 0;
		tube_.back().back().a_concentration = 0.0;
		tube_.back().back().b_concentration = 0.0;
		tube_.back().back().c_concentration = 0.0;
		tube_.back().back().s_concentration = 0.0;


		tube_.back().back().P_block_a = 0.0;
		tube_.back().back().P_block_b = 0.0;
		tube_.back().back().P_block_c = 0.0;



	}


	/*Implements Initial Condition*/
	tube_.back().at(0).a_concentration = a_zero_;
	tube_.back().back().b_concentration = b_zero_;

	/*Suffieciently short times T(!)*/
	tube_.back().at(0).b_concentration = 0;
	tube_.back().back().a_concentration = 0;
	/*----------------*/

	/*Defines Template Vector*/
	tube_template_ = tube_.back();




	return 0;
}
