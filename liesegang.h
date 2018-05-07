/* Source File: Computational Physics 2018
 * Liesegang.h
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#ifndef LIESEGANG_H_
#define LIESEGANG_H_

#include "header.h"




/*Class for Simulation of Liesegang phenomenon*/
class Liesegang {
public:
	Liesegang();
	virtual ~Liesegang();

	/*Tube Unit Space Struct ->
	 * stores value of gas concentrations at each unit of space of the
	 * simulated 'tube'*/
	typedef struct TubeSpaceUnit_{
		double a_concentration;
		double b_concentration;
		double c_concentration;
		double s_concentration;

		/*Keeps track of respective time step;
		 * Note: somewhat ambigous and unnecessary but personally prefered */
		double time_step;

	}TubeSpaceUnit_;


	std::vector<std::vector<TubeSpaceUnit_>> tube;





private:



};

#endif /* LIESEGANG_H_ */
