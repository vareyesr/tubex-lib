/*
 *  CtcExplicitDE
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#ifndef __TUBEX_CTCINTEGROD_H__
#define __TUBEX_CTCINTEGROD_H__

#include "tubex_Ctc.h"
#include "tubex_Slice.h"
#include <vector>
#include <ctime>

namespace tubex
{

	class CtcIntegroD : public Ctc{

	public:

		CtcIntegroD(tubex::Fnc& fnc, Ctc* slice_ctr);

		void contract(TubeVector& x, TubeVector& v,TPropagation t_propa);

	private:

		Ctc* slice_ctr;
		tubex::Fnc& fnc;
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTCINTEGROD_H_ */
