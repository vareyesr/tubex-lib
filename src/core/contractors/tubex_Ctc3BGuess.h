/*
 *  CtcBox class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#ifndef __TUBEX_CTC3BGUESS_H__
#define __TUBEX_CTC3BGUESS_H__

#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include "ibex_Function.h"
#include "tubex_CtcDeriv.h"
#include <vector>
#include <time.h>

namespace tubex
{
	class Ctc3BGuess : public Ctc{
		enum {lb,ub};
		enum {input_gate, output_gate, codomain};
	public:

		Ctc3BGuess(ibex::Fnc& fnc,int bisections=20, double prec=1e-7);

		void contract(TubeVector& x, TubeVector& v, TPropagation t_propa = FORWARD | BACKWARD, int cid_gate=input_gate, bool report=true);

		void ctc_fwdbwd_slices(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos);

		int get_bisections();

		double get_prec();

		void create_slices(Slice & x_slice, std::vector<ibex::Interval> & slices, int cid_gate);

		void change_bisections(int bisections);

		void change_prec(double prec);

		bool oracle_3B(ibex::Interval c1, Slice &xx_slice,Slice& vv_slice,std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int bound, int pos);

		void report(clock_t tStart,TubeVector& x, double old_volume);

	private:
		int bisections;
		double prec;
		ibex::Fnc& fnc;
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTC3BGUESS_H_ */
