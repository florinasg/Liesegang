/* Source File: Computational Physics 2018
 * Liesegang.h
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#ifndef LIESEGANG_H_
#define LIESEGANG_H_

#include "header.h"
#include "Defines.h"




/*Class for Simulation of Liesegang phenomenon*/
class Liesegang {
public:
	Liesegang();
	virtual ~Liesegang();

	/*Tube Unit Space Struct ->
	 * stores value of gas concentrations at each unit of space of the
	 * simulated 'tube'*/

	int InitializeTube();

	int LiesegangGrowth(int mode);

	double Heaviside(double c_conc);





private:

	typedef struct TubeSpaceUnit_{
		TubeSpaceUnit_(): a_concentration(0.0),
				b_concentration(0.0),c_concentration(0.0),
				s_concentration(0.0), time_step(0), position(-1.0), P_block_a(0.0),
				P_block_b(0.0), P_block_c(0.0){};
		double a_concentration;
		double b_concentration;
		double c_concentration;
		double s_concentration;

		/*Keeps track of respective time step;
		 * Note: somewhat ambigous and unnecessary but personally prefered */
		double time_step;

		double position;



		/*----diffusion functional------*/

		/*Introduces Quasi Variables  for Diffusion Functional*/

		double P_block_a;
		double P_block_b;
		double P_block_c;

	}TubeSpaceUnit_;


	/*vector stores vectors of TubeSpaceUnits for every time step
	 * -> length of tube_ then corresponds to the number of timesteps delta_T*/
	std::vector<std::vector<TubeSpaceUnit_>> tube_;

	std::vector<TubeSpaceUnit_> tube_template_; /*Stores Initial Condition & serves as template
	for every succeeding system vector*/

	double delta_X_; /*Self explainatory*/
	double delta_T_; /*Self explainatory*/

	double a_; /*Spatial lower boundary*/
	double b_; /*Spatial upper boundary*/


	double a_zero_;
	double b_zero_;
	double c_zero_;


	double s_zero_;

	double D_a_;
	double D_b_;
	double D_c_;
	double R_;
	double N_one_;
	double N_two_;

	/*Parameters from NON-DIMENSIONALIZATION*/
	double P_ONE_;
	double P_TWO_;


	/*Help Objects*/
	std::ofstream test_file;
	std::ofstream a_concentration_hist;
	std::ofstream b_concentration_hist;
	std::ofstream c_concentration_hist;
	std::ofstream s_concentration_hist;

	double I_time_step_; /*will serve as convenient access index for container*/
	double alpha_a_;
	double alpha_b_;
	double alpha_c_;

};

#endif /* LIESEGANG_H_ */
