#include "tubex_CtcDynBasic.h"

using namespace std;
using namespace ibex;


namespace tubex
{
	CtcDynBasic::CtcDynBasic(ibex::Fnc& fnc, double prec): fnc(fnc), prec(prec)
	{
		/*check input*/
		assert(prec >= 0);
	}

	void CtcDynBasic::contract(Slice& x, Slice& v, TPropagation t_propa)
	{
		/*check if the domain of x is the same as v*/
		assert(x.domain() == v.domain());
//		assert(TubeVector::same_slicing(x, v));

		if(m_fast_mode)
			ctc_deriv.set_fast_mode(true);

		ctc_deriv.contract(x,v,t_propa);
		ctc_fwd(x,v);
	}


	void CtcDynBasic::ctc_fwd(Slice &x, Slice &v)
	{
		/*envelope*/
//		v.set_envelope(fnc.eval_vector(x.codomain)[pos]);
	}

	double CtcDynBasic::get_prec()
	{
		return this->prec;
	}

	void CtcDynBasic::set_prec(double prec)
	{
		this->prec = prec;
	}
}



