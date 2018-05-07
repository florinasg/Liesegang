/* Source File: Computational Physics 2018
 * header.h
 *
 *  Created on: 07.05.2018
 *      Author: Sara Gasparini
 */

#ifndef HEADER_H_
#define HEADER_H_


#include <fstream>
#include <armadillo>
#include <math.h>
#include <cmath>
#include <iostream>
#include <vector>



/*easy console output*/
template <typename T>
void print(T t)
{

  std::cout << t << " " ;
}

template<typename T, typename... Args>
void print(T t, Args... args)
{

  std::cout << t << " ";
  print(args...) ;
}




#endif /* HEADER_H_ */
