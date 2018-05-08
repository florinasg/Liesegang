/* Source File: Computational Physics 2018
 * main.cpp
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#include "liesegang.h"
#include "Defines.h"


/*args = [
 * -
 * -
 * -
 *
 * ]*/
int main(int args, char * argv[])
{

	/*EASY: Create Instance wit Defines
	 * -> ADVANCED: Use argv[]*/
	Liesegang *liesegang = new Liesegang();
	liesegang->InitializeTube();
	liesegang->LiesegangGrowth(true);

	return 0;
}


