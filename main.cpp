/* Source File: Computational Physics 2018
 * main.cpp
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#include "liesegang.h"
#include "Defines.h"
#include <istream>
#include <iostream>


/*args = [
 * -
 * - NO ARGS
 * -
 *
 * ]*/
int main(int args, char * argv[])
{

	/*EASY: Create Instance wit Defines
	 * -> ADVANCED: Use argv[]*/
	Liesegang *liesegang = new Liesegang();
	liesegang->InitializeTube();





	if(args == 1)
	{
		std::cout << "Liesegang Growth with Constant Diffusion Constant..." << std::endl;
		liesegang->LiesegangGrowth(0);
		std::cout << "Admire the Structure..";
	}
	else if(atoi(argv[1]) == 0 || atoi(argv[1]) == 1)
	{
		std::cout << "Liesegang Growth with Diffusion Functionsl..." << std::endl;
		liesegang->LiesegangGrowth(atoi(argv[1]));
		std::cout << "Admire the Structure..";
	}

	else
	{
		std::cout << "ERROR: Please call exe either with [0,1] or no argument" << std::endl;
	}



	return 0;
}


