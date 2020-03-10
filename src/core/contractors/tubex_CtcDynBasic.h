/*
 * tubex_CtcDynBasic.h
 *
 *  Created on: Mar 6, 2020
 *      Author: Victor Reyes
 *      		Gilles Trombettoni
 */

#ifndef __TUBEX_CTCDYNBASIC_H__
#define __TUBEX_CTCDYNBASIC_H__

#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include "ibex_Function.h"
#include "tubex_CtcDeriv.h"
#include <vector>
#include <time.h>

namespace tubex
{

	class CtcDynBasic : public Ctc{

	public:

		CtcDynBasic(ibex::Fnc& fnc, double prec = 0.);
		/*
		 * This method performs a contraction for the TubeVector x.
		 * Note that the timesteps between the Tubes of x must be identically the same.
		 */
		void contract(Slice& x, Slice& v, TPropagation t_propa);
		/*
		 * ctc_fwd manages to make an evaluation of the current Slice in order to contract and update v
		 */
		void ctc_fwd(Slice &x, Slice &v);
		/*
		 *  used to obtain the current precision of the iterative method.
		 */
		double get_prec();
		/*
		 * changes the value of the precision
		 */
		void set_prec(double prec);

	private:
		CtcDeriv ctc_deriv;
		ibex::Fnc& fnc;
		double prec;
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTCDYNBASIC_H_ */

