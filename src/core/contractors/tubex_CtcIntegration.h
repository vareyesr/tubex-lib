/*
 *  CtcIntegration class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#ifndef __TUBEX_CTCINTEGRATION_H__
#define __TUBEX_CTCINTEGRATION_H__

#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include "ibex_Function.h"
#include "tubex_CtcDeriv.h"
#include <vector>
#include <ctime>


namespace tubex
{

	class CtcIntegration : public Ctc{

	public:

		/*
		 * CtcIntegration is a contractor that works at the Tube level. It requires as input an evolution function and a Slice contractor.
		 * Currently, the slice contractors that can be used are: CtcDynCid and CtcDynCidGuess
		 */
		CtcIntegration(ibex::Fnc& fnc, Ctc* slice_ctr);
		/*
		 * This method performs a contraction for the TubeVector x.
		 * Note that the timesteps between the Tubes of x must be identically the same.
		 */
		void contract(TubeVector& x, TubeVector& v, double time_dom, TPropagation t_propa, bool report=false);
		/*
		 * Same as above, but v is computed inside.
		 */
		void contract(TubeVector& x, double time_dom, TPropagation t_propa, bool report=false);
		/*
		 * ctc_fwd manages to make an evaluation of the current Slice in order to contract and update v
		 */
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos);

		/*
		 * used to report the results
		 */
		void report(clock_t tStart,TubeVector& x, double old_volume);

		/*
		 * returns the final time reached during contraction
		 */
		double get_finaltime();

	private:
		Ctc* slice_ctr;
		ibex::Fnc& fnc;
		double finaltime;
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTCINTEGRATION_H_ */
