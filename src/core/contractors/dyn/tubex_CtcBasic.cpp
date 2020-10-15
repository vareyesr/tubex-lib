/*
*  Author(s):    Victor Reyes, Gilles Trombettoni
*  Created:      Mar 6, 2020
*  Updated:      Sep 11, 2020
*/

#include "tubex_CtcBasic.h"

using namespace std;
using namespace ibex;


namespace tubex
{
	CtcBasic::CtcBasic(const TFnc& fnc, double prec): fnc(fnc), prec(prec)
	{
		/*check input*/
		assert(prec >= 0);
	}

	bool CtcBasic::contract(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa)
	{
		/*check if the domains are the same*/
		Interval to_try(x_slice[0]->tdomain());
		for (auto xslice : x_slice)
			assert(to_try == xslice->tdomain());

		//check if the gates used to contract are bounded
		for (auto xslice : x_slice){
			if ((t_propa & TimePropag::FORWARD) && (xslice->input_gate().is_unbounded()))
				return false;
			else if ((t_propa & TimePropag::BACKWARD) && (xslice->output_gate().is_unbounded()))
				return false;
		}

		bool fix_point_n;
		bool first_iteration = true;

	do{
		fix_point_n = false;
		double volume = 0.;

		/*volume before contraction*/
		for (auto xslice : x_slice)
		  volume+=xslice->volume();

		for (unsigned int i = 0 ; i < x_slice.size() ; i++){

			Slice aux_slice_x(*x_slice[i]);
			Slice aux_slice_v(*v_slice[i]);

			double sx;
			/*without polygons*/
			if(m_fast_mode)
				ctc_deriv.set_fast_mode(true);
			if (get_reasoning_slice()){
				do
				{
					sx = aux_slice_x.volume();
					ctc_deriv.contract(aux_slice_x, aux_slice_v,t_propa);
					ctc_fwd(aux_slice_x, aux_slice_v, x_slice, v_slice, i);
				} while ((1-(aux_slice_x.volume()/sx)) > get_prec());
			}
			else{
				ctc_deriv.contract(aux_slice_x, aux_slice_v,t_propa);
				ctc_fwd(aux_slice_x, aux_slice_v, x_slice, v_slice, i);
			}

			/*Replacing the old domains with the new ones*/
			x_slice[i]->set_envelope(aux_slice_x.codomain()); v_slice[i]->set_envelope(aux_slice_v.codomain());
			x_slice[i]->set_input_gate(aux_slice_x.input_gate()); v_slice[i]->set_input_gate(aux_slice_v.input_gate());
			x_slice[i]->set_output_gate(aux_slice_x.output_gate()); v_slice[i]->set_output_gate(aux_slice_v.output_gate());
		}

		/*volume after contraction*/
		double volume_after=0.;
		for (auto xslice : x_slice)
		  volume_after+=xslice->volume();


		if (volume_after < volume)
		  fix_point_n = true;

		if ((first_iteration) && !(fix_point_n))
		  return false;

		first_iteration = false;

		if (!get_reasoning_slice()) fix_point_n=false;

	} while(fix_point_n);

		return true;
	}

	void CtcBasic::contract(std::vector<Domain*>& v_domains){
		//todo is it useful to implement it?
	}

	void CtcBasic::ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, unsigned int pos)
	{
		/*envelope*/
		IntervalVector envelope(x_slice.size()+1);
		envelope[0] = x.tdomain();

		for (unsigned int i = 0 ; i < x_slice.size() ; i++){
			if (i==pos)
				envelope[i+1] = x.codomain();
			else
				envelope[i+1] = x_slice[i]->codomain();
		}
		v.set_envelope(fnc.eval_vector(envelope)[pos]);
	}

	double CtcBasic::get_prec()
	{
		return this->prec;
	}

	void CtcBasic::set_prec(double prec)
	{
		this->prec = prec;
	}

	bool CtcBasic::get_reasoning_slice(){
		return this->m_reasoning_slice;
	}

	void CtcBasic::set_reasoning_slice(bool reasoning_slice){
		this->m_reasoning_slice = reasoning_slice;
	}
}
