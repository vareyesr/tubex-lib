/*
 *  Author(s):    Victor Reyes, Gilles Trombettoni
 *  Created:      Mar 6, 2020
 *  Updated:      Sep 11, 2020
 */

#include "tubex_CtcCid.h"
#include <deque>


using namespace std;
using namespace ibex;


namespace tubex
{
	CtcCid::CtcCid(const TFnc& fnc,int scid, double prec): fnc(fnc), scid(scid), prec(prec)
	{
		/*check inputs*/
		assert(scid > 0.);
		assert(prec >= 0);
	}

	bool CtcCid::contract(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa)
	{
		//checks that the domain of each slice is the same.
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
		//
		//polygone method setter
		if(m_fast_mode)
			ctc_deriv.set_fast_mode(true);

		//save volume before contraction
		double volume_old = 0;
		for (auto xslice : x_slice)
			volume_old = volume_old + xslice->volume();

		/*propagation engine*/
		bool fix_point;
		do{
			fix_point = false;
			double volume_1 = 0; double volume_2 = 0;
			for (auto xslice : x_slice)
				volume_1 = volume_1 + xslice->volume();

			if (get_propagation_engine() == 0)
				AtomicPropagationEngine(x_slice,v_slice,t_propa);
			else if (get_propagation_engine() == 1)
				FullPropagationEngine(x_slice,v_slice,t_propa);
			for (auto xslice : x_slice)
				volume_2= volume_2 + xslice->volume();
			if (1-(volume_2/volume_1) > get_prec()) fix_point = true;
		} while (fix_point);

		/*empty test*/
		for (auto xslice : x_slice)
			if (xslice->is_empty()) return false;

		/*incrementality test*/
		double volume_new = 0;
		for (auto xslice : x_slice)
			volume_new = volume_new + xslice->volume();

		if (volume_old == volume_new) return false;

		return true;
	}

	void CtcCid::ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos)
	{
		/*envelope*/
		IntervalVector envelope(x_slice.size()+1);
		envelope[0] = x.tdomain();

		for (unsigned int i = 0 ; i < x_slice.size() ; i++){
			if (i==(unsigned)pos)
				envelope[i+1] = x.codomain();
			else
				envelope[i+1] = x_slice[i]->codomain();
		}
		v.set_envelope(fnc.eval_vector(envelope)[pos]);
	}


	void CtcCid::ctc_fwd(Slice &x, Slice &v, std::vector<Slice> x_slice, std::vector<Slice> v_slice, int pos)
	{
		/*envelope*/
		IntervalVector envelope(x_slice.size()+1);
		envelope[0] = x.tdomain();

		for (unsigned int i = 0 ; i < x_slice.size() ; i++){
			if (i==(unsigned)pos)
				envelope[i+1] = x.codomain();
			else
				envelope[i+1] = x_slice[i].codomain();
		}

		v.set_envelope(fnc.eval_vector(envelope)[pos]);
	}

	double CtcCid::get_scid()
	{
		return this->scid;
	}

	double CtcCid::get_prec()
	{
		return this->prec;
	}

	int CtcCid::get_propagation_engine(){
		return this->engine;
	}

	void CtcCid::set_scid(int scid)
	{
		this->scid = scid;
	}

	void CtcCid::set_prec(double prec)
	{
		this->prec = prec;
	}

	void CtcCid::set_propagation_engine(int engine){
		this->engine = engine;
		if (this->engine == 0) this->set_prec(0);      //difficult to check, as the contraction is more volatile.
		else if (this->engine == 1) this->set_prec(0.01); // todo: seems a good tradeoff between cpu-time and contraction?
	}



	void CtcCid::FullPropagationEngine(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa){

		/*create the contractor queue: format contraint - variable, 0: for ctc_deriv, 1 for fwd*/
		std::deque< vector<int> > contractorQ;

		vector<int> ctr_var;
		for (int i = 0 ; i < 2 ; i++) ctr_var.push_back(-1);

		/*Initialization contractor array - bool*/
		vector<bool > isPresent;

		for (unsigned int i = 0 ; i < x_slice.size() ; i++)
			isPresent.push_back(false);


		/*going throw all the variables*/
		for (unsigned int i = 0 ; i < x_slice.size() ; i++){

			std::vector<ibex::Interval> x_subslices;
			x_subslices.clear();
			/*create the sub-slices*/
			create_subintervals(*x_slice[i],x_subslices, t_propa);

			/*Hull for each dimension on x and v*/
			vector<Interval> hull_input_x; vector<Interval> hull_input_v;
			vector<Interval> hull_output_x; vector<Interval> hull_output_v;
			vector<Interval> hull_codomain_x; vector<Interval> hull_codomain_v;


			/*Initiliazation*/
			for (unsigned int j = 0 ; j < x_slice.size() ; j++){
				hull_input_x.push_back(Interval::EMPTY_SET); hull_input_v.push_back(Interval::EMPTY_SET);
				hull_output_x.push_back(Interval::EMPTY_SET); hull_output_v.push_back(Interval::EMPTY_SET);
				hull_codomain_x.push_back(Interval::EMPTY_SET); hull_codomain_v.push_back(Interval::EMPTY_SET);
			}

			for (unsigned int k = 0 ; k < x_subslices.size() ; k++){

				/*restore with the current domains*/
				vector<Slice> aux_slice_x; aux_slice_x.clear();
				vector<Slice> aux_slice_v; aux_slice_v.clear();

				for (unsigned int j = 0 ; j < x_slice.size() ; j++){
					aux_slice_x.push_back(*x_slice[j]);
					aux_slice_v.push_back(*v_slice[j]);
				}

				/*Set the gate depending on the direction of the contraction*/
				if (t_propa & TimePropag::FORWARD)
					aux_slice_x[i].set_input_gate(x_subslices[k]);
				else if (t_propa & TimePropag::BACKWARD)
					aux_slice_x[i].set_output_gate(x_subslices[k]);

				/*push the first element to the contractor queue*/
				ctr_var[0] = 0; ctr_var[1] = i;
				contractorQ.push_back(ctr_var);

				/*FIFO queue*/
				do{
					/*get what contractor should be called*/
					int contractor = contractorQ.front()[0];
					/*get the variable that is going to be contracted*/
					int variable = contractorQ.front()[1];
					isPresent[variable] = false;
					/*pop the first element*/
					contractorQ.pop_front();
					/*contract*/
					if (contractor == 0 ){ //call ctc_deriv
						/*save the corresponding domain*/
						double size_x = aux_slice_x[variable].codomain().diam();
						ctc_deriv.contract(aux_slice_x[variable],aux_slice_v[variable],t_propa);

						if ((1-(aux_slice_x[variable].codomain().diam()/size_x)) > this->get_prec()){
							ctr_var[0] = 1; ctr_var[1] = variable;
							contractorQ.push_front(ctr_var);
						}
					}
					else if (contractor == 1){ //call ctc_fwd
						/*save the corresponding domain*/
						double size_v = aux_slice_v[variable].codomain().diam();
						ctc_fwd(aux_slice_x[variable], aux_slice_v[variable], aux_slice_x, aux_slice_v, variable);

						/*add the contraints not included in isPresent*/
						if ((1-(aux_slice_v[variable].codomain().diam()/size_v)) > this->get_prec()){
							for (unsigned int j = 0 ; j < x_slice.size() ; j++ ){
								if ((!isPresent[j]) && (j!=(unsigned)variable)){
									ctr_var[0] = 1; ctr_var[1] = j;
									contractorQ.push_front(ctr_var);
									isPresent[j] = true;
								}
							}
							ctr_var[0] = 0; ctr_var[1] = variable;
							contractorQ.push_front(ctr_var);
						}
					}

				} while (contractorQ.size() > 0);

				/*The union of the current Slice is made*/
				for (unsigned int j = 0 ; j < x_slice.size() ; j++){
					hull_input_x[j] |= aux_slice_x[j].input_gate(); hull_input_v[j] |= aux_slice_v[j].input_gate();
					hull_output_x[j] |= aux_slice_x[j].output_gate(); hull_output_v[j] |= aux_slice_v[j].output_gate();
					hull_codomain_x[j] |= aux_slice_x[j].codomain(); hull_codomain_v[j] |= aux_slice_v[j].codomain();
				}
			}

			/*replacing the old domains*/
			for (unsigned int j = 0 ; j < x_slice.size() ; j++){
				x_slice[j]->set_envelope(x_slice[j]->codomain() & hull_codomain_x[j]); v_slice[j]->set_envelope(v_slice[j]->codomain() & hull_codomain_v[j]);
				x_slice[j]->set_input_gate(x_slice[j]->input_gate() & hull_input_x[j]); v_slice[j]->set_input_gate(v_slice[j]->input_gate() & hull_input_v[j]);
				x_slice[j]->set_output_gate(x_slice[j]->output_gate() & hull_output_x[j]); v_slice[j]->set_output_gate(v_slice[j]->output_gate() & hull_output_v[j]);
			}
		}
	}

	void CtcCid::AtomicPropagationEngine(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa){

		for (unsigned int i = 0 ; i < x_slice.size() ; i++){

			std::vector<ibex::Interval> x_subslices;
			x_subslices.clear();

			/*create the sub-slices*/
			create_subintervals(*x_slice[i],x_subslices, t_propa);

			/*For each slice on $t$ compute the corresponding the hull */
			Interval hull_input_x = Interval::EMPTY_SET; Interval hull_input_v = Interval::EMPTY_SET;
			Interval hull_output_x = Interval::EMPTY_SET; Interval hull_output_v = Interval::EMPTY_SET;
			Interval hull_codomain_x = Interval::EMPTY_SET; Interval hull_codomain_v = Interval::EMPTY_SET;

			/*treat each subslice, make the hull and then intersect*/
			for (unsigned int j = 0 ; j < x_subslices.size() ; j++){

				/*Temporal slices on $x$ and $v$*/
				Slice aux_slice_x(*x_slice[i]);
				Slice aux_slice_v(*v_slice[i]);

				if (t_propa & TimePropag::FORWARD)
					aux_slice_x.set_input_gate(x_subslices[j]);
				else if (t_propa & TimePropag::BACKWARD)
					aux_slice_x.set_output_gate(x_subslices[j]);

				/*Fixpoint for each sub-slice at each tube*/
				double sx;

				/*without polygons*/
				if(m_fast_mode)
					ctc_deriv.set_fast_mode(true);

				do
				{
					sx = aux_slice_x.volume();
					ctc_deriv.contract(aux_slice_x, aux_slice_v,t_propa);
					ctc_fwd(aux_slice_x, aux_slice_v, x_slice, v_slice, i);

				} while((1-(aux_slice_x.volume()/sx)) > get_prec());

				/*The union of the current Slice is made.*/
				hull_input_x |= aux_slice_x.input_gate(); hull_input_v |= aux_slice_v.input_gate();
				hull_output_x |= aux_slice_x.output_gate(); hull_output_v |= aux_slice_v.output_gate();
				hull_codomain_x |= aux_slice_x.codomain(); hull_codomain_v |= aux_slice_v.codomain();
			}

			/*Replacing the old domains with the new ones*/

			x_slice[i]->set_envelope(hull_codomain_x); v_slice[i]->set_envelope(hull_codomain_v);
			x_slice[i]->set_input_gate(hull_input_x); v_slice[i]->set_input_gate(hull_input_v);
			x_slice[i]->set_output_gate(hull_output_x); v_slice[i]->set_output_gate(hull_output_v);
		}
	}

	void CtcCid::create_subintervals(Slice& x_slice, std::vector<ibex::Interval> & x_slices, TimePropag t_propa)
	{
		/*Varcid in the input gate*/
		if (t_propa & TimePropag::FORWARD){
			double size_interval = x_slice.input_gate().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.input_gate().lb()+i*size_interval,x_slice.input_gate().lb()+size_interval*(i+1)));
			}
		}

		/*Varcid in the output gate*/
		else if (t_propa & TimePropag::BACKWARD){
			double size_interval = x_slice.output_gate().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.output_gate().lb()+i*size_interval,x_slice.output_gate().lb()+size_interval*(i+1)));
			}
		}
	}
}

