/* Source File: Computational Physics 2018
 * heaviside.cpp
 *
 *  Created on: 08.05.2018
 *      Author: Sara Gasparini
 */


#include "liesegang.h"


double Liesegang::Heaviside(double c_conc)
{
 if(c_conc < c_zero_)
 {
	 return 0;
 }

 else if (c_conc == c_zero_)
 {
	 return 1/2;
 }

 else if(c_conc >  c_zero_)
 {
	 return 1;
 }

return 0.0;
}

