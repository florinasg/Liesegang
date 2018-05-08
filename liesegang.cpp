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
	D_c_= double(DEF_D_C);
	R_= double(DEF_R);
	N_one_ = double(DEF_N);
	N_two_= double(DEF_N) ;


	/*TODO: For Later Use*/
	P_ONE_= double() ;
	P_TWO_= double() ;

	alpha_a_ = alpha_b_ = alpha_c_ = 0.0;

	delta_T_ = 0.0;
	delta_X_ = 0.0;


	std::string prefix = "(Da="+std::to_string(D_a_)+",Db="+std::to_string(D_b_)
	+",Dc="+ std::to_string(D_c_)+",N1="+std::to_string(N_one_)+",N2="+std::to_string(N_two_)+
	",R="+std::to_string(R_)+",a0="+std::to_string(a_zero_)+",b0="+
	std::to_string(b_zero_)+",c0="+std::to_string(c_zero_)+")";

	/*Creates Directory for */
	std::string dir_name =".\\Results_"+prefix;
	int check = _mkdir(dir_name.c_str());

	a_concentration_hist.open("A_Concentration_"+prefix+".csv",std::fstream::trunc);
	b_concentration_hist.open("B_Concentration_"+prefix+".csv",std::fstream::trunc);
	c_concentration_hist.open("C_Concentration_"+prefix+".csv",std::fstream::trunc);
	s_concentration_hist.open("S_Concentration_"+prefix+".csv",std::fstream::trunc);







}

Liesegang::~Liesegang()
{
	a_concentration_hist.close();
	b_concentration_hist.close();
	c_concentration_hist.close();
	s_concentration_hist.close();

	test_file.close();
}



