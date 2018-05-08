/* Source File: Computational Physics 2018
 * Liesegang.cpp
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#include "Liesegang.h"

Liesegang::Liesegang() {


	test_file.open("Control_File.csv");

	/*Initialization of Parameter*/

	a_ = double(DEF_LOWER_BOUND) ;
	b_= double(DEF_UPPER_BOUND) ;


	a_zero_= double(DEF_A_ZERO) ;
	b_zero_= double(DEF_B_ZERO) ;
	c_zero_= double(DEF_C_ZERO) ;

	D_a_= double(DEF_D_A) ;
	D_b_= double(DEF_D_B) ;
	D_c_= double(DEF_C_ZERO) ;
	R_= double(DEF_R);
	N_one_ = double(DEF_N);
	N_two_= double(DEF_N) ;


	/*TODO: For Later Use*/
	P_ONE_= double() ;
	P_TWO_= double() ;



}

Liesegang::~Liesegang()
{
	test_file.close();
}



