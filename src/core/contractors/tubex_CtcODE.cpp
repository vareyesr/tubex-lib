/*
 *  CtcDynCid Class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#include "tubex_CtcODE.h"

using namespace std;
using namespace ibex;


namespace tubex
{
	CtcODE::CtcODE(std::vector<double> observations, CtcIntegration integration): observations(observations), integration(integration)
	{
		/*check input, size must be greater than 0*/
		assert(observations.size() > 0.);
		/*check if it is sorted?*/

	}

	void CtcODE::contract(TubeVector& x, TubeVector& v)
	{

		/*FORWARD Phase*/
		for (int i = 0 ; i < observations.size() ; i++){
			integration.contract(x,v,observations[i],FORWARD);
			int ii = 0;
			for (int j = i+1 ; j < observations.size() ; j++){
				if (integration.get_finaltime() > observations[j])
					ii++;
			}
			i+=ii;
		}

		/*BACKWARD Phase*/
		for (int i = observations.size()-1 ; i >= 0 ; i--){
			integration.contract(x,v,observations[i],BACKWARD);
			int ii = 0;
			for (int j = i ; j >=0 ; j--){
				if (integration.get_finaltime() < observations[j])
					ii--;
			}
			i+=ii;
		}
	}
}
