/*
 *  CtcBox class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes
 */
#include "tubex_CtcDynCidGuess.h"

using namespace std;
using namespace ibex;



namespace tubex
{
	CtcDynCidGuess::CtcDynCidGuess(ibex::Fnc& fnc,int bisections, double prec): fnc(fnc), bisections(bisections), prec(prec)
	{
		assert(bisections >= 0.);
		assert(prec >= 0);
	}


	bool CtcDynCidGuess::contract(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TPropagation t_propa)
	{
		//checks that the domain of each slice is the same.
		Interval to_try(x_slice[0]->domain());
		double total_volume = 0;

		for (int i = 1 ; i < x_slice.size(); i++){
			assert(to_try == x_slice[i]->domain());
			total_volume+=x_slice[i]->volume();
		}

		bool fix_point_l;
		bool first_iteration = true;
		do{
			fix_point_l = false;
			/*with the bounds of each variable compute the impact on the output door*/
			vector<Slice> x_slice_bounds;
			vector<Slice> v_slice_bounds;
			x_slice_bounds.clear(); v_slice_bounds.clear();

			for (int i = 0 ; i < x_slice.size() ; i++){
				x_slice_bounds.push_back(*x_slice[i]);
				v_slice_bounds.push_back(*v_slice[i]);
			}

			/*Guess procedure, just using the bounds*/
			bool fix_point_n;
			do{
				fix_point_n = false;
				/*work for each dimension*/
				for (int i = 0 ; i < x_slice.size() ;i++){

					std::vector<ibex::Interval> x_subslices;
					x_subslices.clear();
					create_slices(x_slice_bounds[i],x_subslices, t_propa);

					Interval hull_input_x = Interval::EMPTY_SET; Interval hull_input_v = Interval::EMPTY_SET;
					Interval hull_output_x = Interval::EMPTY_SET; Interval hull_output_v = Interval::EMPTY_SET;
					Interval hull_codomain_x = Interval::EMPTY_SET; Interval hull_codomain_v = Interval::EMPTY_SET;

					for (int j = 0 ; j < x_subslices.size() ; j++){

						/*Temporal slices on $x$ and $v$*/
						Slice aux_slice_x(x_slice_bounds[i]);
						Slice aux_slice_v(v_slice_bounds[i]);

						if (t_propa & FORWARD)
							aux_slice_x.set_input_gate(x_subslices[j]);
						else if (t_propa & BACKWARD)
							aux_slice_x.set_output_gate(x_subslices[j]);

						/*Fixpoint for each sub-slice at each tube*/
						double sx;

						/*without polygons*/
						if (m_fast_mode)
							ctc_deriv.set_fast_mode(true);
						do
						{
							sx = aux_slice_x.volume();
							ctc_deriv.contract(aux_slice_x, aux_slice_v,t_propa);
							ctc_fwd(aux_slice_x, aux_slice_v, x_slice, v_slice, i);

						} while(sx-aux_slice_x.volume()>get_prec());


						/*The union of the current guess is made.*/
						hull_input_x |= aux_slice_x.input_gate(); hull_input_v |= aux_slice_v.input_gate();
						hull_output_x |= aux_slice_x.output_gate(); hull_output_v |= aux_slice_v.output_gate();
						hull_codomain_x |= aux_slice_x.codomain(); hull_codomain_v |= aux_slice_v.codomain();

					}

					double volume = x_slice_bounds[i].volume();

					/*Replacing the slice with impact of the output (either the input) gates in all the dimensions*/
					x_slice_bounds[i].set_envelope(hull_codomain_x);  v_slice_bounds[i].set_envelope(hull_codomain_v);
					x_slice_bounds[i].set_input_gate(hull_input_x); v_slice_bounds[i].set_input_gate(hull_input_v);
					x_slice_bounds[i].set_output_gate(hull_output_x); v_slice_bounds[i].set_output_gate(hull_output_v);


					if (volume > x_slice_bounds[i].volume())
						fix_point_n = true;

				}
			} while(fix_point_n);

			/*3B part*/
			for (int i = 0 ; i < x_slice.size() ; i++){
				double aux_envelope = x_slice[i]->codomain().diam();
				Interval remove_ub,remove_lb;
				if (t_propa & FORWARD){
					remove_ub = Interval(x_slice_bounds[i].output_gate().ub(),x_slice[i]->output_gate().ub());
					remove_lb = Interval(x_slice[i]->output_gate().lb(),x_slice_bounds[i].output_gate().lb());
				}
				else if (t_propa & BACKWARD){
					remove_ub = Interval(x_slice_bounds[i].input_gate().ub(),x_slice[i]->input_gate().ub());
					remove_lb = Interval(x_slice[i]->input_gate().lb(),x_slice_bounds[i].input_gate().lb());
				}

				//todo: check if necessary
				if (remove_ub.diam()>0){
					var3Bcheck(remove_ub,ub,i,x_slice,v_slice,t_propa);
					ctc_deriv.contract(*x_slice[i], *v_slice[i],t_propa);
					ctc_fwd(*x_slice[i], *v_slice[i], x_slice, v_slice, i);
					ctc_deriv.contract(*x_slice[i], *v_slice[i],t_propa);
				}
				if (remove_lb.diam()>0){
					var3Bcheck(remove_lb,lb,i,x_slice,v_slice,t_propa);
					ctc_deriv.contract(*x_slice[i], *v_slice[i],t_propa);
					ctc_fwd(*x_slice[i], *v_slice[i], x_slice, v_slice, i);
					ctc_deriv.contract(*x_slice[i], *v_slice[i],t_propa);
				}

				/*empty test*/
				if (x_slice[i]->is_empty()) return false;
				if (aux_envelope > x_slice[i]->codomain().diam())
					fix_point_l = true;
			}
			if ((first_iteration) && !(fix_point_l))
					return false;

			first_iteration = false;
		} while (fix_point_l);

		return true;
	}

	void CtcDynCidGuess::ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos)
	{
		/*envelope*/
		IntervalVector envelope(x_slice.size());

		for (int i = 0 ; i < x_slice.size() ; i++){
			if (i==pos)
				envelope[i] = x.codomain();
			else
				envelope[i] = x_slice[i]->codomain();
		}
			v.set_envelope(fnc.eval_vector(envelope)[pos]);
	}

	int CtcDynCidGuess::get_bisections()
	{
		return this->bisections;
	}

	double CtcDynCidGuess::get_prec()
	{
		return this->prec;
	}

	void CtcDynCidGuess::change_prec(double prec){
		this->prec = prec;
	}

	void CtcDynCidGuess::change_bisections(int bisections)
	{
		this->bisections = bisections;
	}

	void CtcDynCidGuess::create_slices(Slice& x_slice, std::vector<ibex::Interval> & x_slices, TPropagation t_propa)
	{

		/*Varcid in the input gate*/
		if (t_propa & FORWARD){
			if (x_slice.input_gate().diam() != 0){
				x_slices.push_back(Interval(x_slice.input_gate().lb()));
				x_slices.push_back(Interval(x_slice.input_gate().ub()));
			}
		}

		/*Varcid in the output gate*/
		else if (t_propa & BACKWARD){
			if (x_slice.output_gate().diam() == 0)
				x_slices.push_back(Interval(x_slice.output_gate().ub()));
			else{
				x_slices.push_back(Interval(x_slice.output_gate().lb()));
				x_slices.push_back(Interval(x_slice.output_gate().ub()));
			}
		}
	}

	void CtcDynCidGuess::var3Bcheck(ibex::Interval remove_bound ,int bound, int pos ,std::vector<Slice*> & x_slice, std::vector<Slice*> v_slice, TPropagation t_propa)
	{

		ctc_deriv.set_fast_mode(false);
		bool fix_point_n;
		vector<Slice> x_slice_bounds;
		vector<Slice> v_slice_bounds;
		for (int i = 0 ; i < x_slice.size() ; i++){
			x_slice_bounds.push_back(*x_slice[i]);
			v_slice_bounds.push_back(*v_slice[i]);
		}
		if (t_propa & FORWARD)
			x_slice_bounds[pos].set_output_gate(remove_bound);
		else if (t_propa & BACKWARD)
			x_slice_bounds[pos].set_input_gate(remove_bound);
		/*try to remove the complete interval*/
		for (int i = 0 ;  i < x_slice.size() ; i++){
			double sx;

			do
			{
				sx = x_slice_bounds[i].volume();
				ctc_deriv.contract(x_slice_bounds[i], v_slice_bounds[i],t_propa);
				ctc_fwd(x_slice_bounds[i], v_slice_bounds[i], x_slice, v_slice, i);
			} while(sx-x_slice_bounds[i].volume()>get_prec());

			/*if something is empty means that we can remove the complete interval*/
			if (x_slice_bounds[i].is_empty()){
				if (t_propa & FORWARD){
					if (bound == ub)
						x_slice[pos]->set_output_gate(Interval(x_slice[pos]->output_gate().lb(),remove_bound.lb()));
					else if (bound == lb)
						x_slice[pos]->set_output_gate(Interval(remove_bound.ub(),x_slice[pos]->output_gate().ub()));
				}
				if (t_propa & BACKWARD){
					if (bound == ub)
						x_slice[pos]->set_input_gate(Interval(x_slice[pos]->input_gate().lb(),remove_bound.lb()));
					else if (bound == lb)
						x_slice[pos]->set_input_gate(Interval(remove_bound.ub(),x_slice[pos]->input_gate().ub()));
				}
				return;
			}
		}


		/*dichotomic way*/
		Interval half;bool success;
//		if (t_propa & FORWARD){
//			if (bound == ub)
//				remove_bound = Interval(remove_bound.lb(),x_slice_bounds[pos].output_gate().ub());
//			else if (bound == lb)
//				remove_bound = Interval(x_slice_bounds[pos].output_gate().lb(),remove_bound.ub());
//		}
//		else if (t_propa & BACKWARD){
//			if (bound == ub)
//				remove_bound = Interval(remove_bound.lb(),x_slice_bounds[pos].input_gate().ub());
//			else if (bound == lb)
//				remove_bound = Interval(x_slice_bounds[pos].input_gate().lb(),remove_bound.ub());
//		}
		for (int k = 0 ;  k < 10 ; k++){
			/*restore domains*/
			x_slice_bounds.clear(); v_slice_bounds.clear();
			for (int i = 0 ; i < x_slice.size() ; i++){
				x_slice_bounds.push_back(*x_slice[i]);
				v_slice_bounds.push_back(*v_slice[i]);
			}

			if (bound == ub)
				half = Interval(remove_bound.mid(),remove_bound.ub());

			else if (bound == lb)
				half = Interval(remove_bound.lb(),remove_bound.mid());
			if (t_propa & FORWARD)
				x_slice_bounds[pos].set_output_gate(half);
			else if (t_propa & BACKWARD)
				x_slice_bounds[pos].set_input_gate(half);
			success = false;
			for (int i = 0 ;  i < x_slice_bounds.size() ; i++){
				double sx;
				CtcDeriv ctc_deriv;
//				/*without polygons*/
				do
				{
					sx = x_slice_bounds[i].volume();
					ctc_deriv.contract(x_slice_bounds[i], v_slice_bounds[i],t_propa);
					ctc_fwd(x_slice_bounds[i], v_slice_bounds[i], x_slice, v_slice, i);
				} while(sx-x_slice_bounds[i].volume()>get_prec());
				/*if something is empty means that we can remove the half*/

				if (x_slice_bounds[i].is_empty()){
					success = true;
					if (bound == ub)
						remove_bound = Interval(remove_bound.lb(),remove_bound.mid());
					else if (bound == lb)
						remove_bound = Interval(remove_bound.mid(),remove_bound.ub());
				}
				if (success)
					break;
			}
			/*cant remove half, finish!*/
			if (!success){
				break;
			}

		}

		if (t_propa & FORWARD){
			if (bound==ub)
				x_slice[pos]->set_output_gate(Interval(x_slice[pos]->output_gate().lb(),remove_bound.ub()));
			else if (bound==lb)
				x_slice[pos]->set_output_gate(Interval(remove_bound.lb(),x_slice[pos]->output_gate().ub()));
		}
		else if (t_propa & BACKWARD){
			if (bound==ub)
				x_slice[pos]->set_input_gate(Interval(x_slice[pos]->input_gate().lb(),remove_bound.ub()));
			else if (bound==lb)
				x_slice[pos]->set_input_gate(Interval(remove_bound.lb(),x_slice[pos]->input_gate().ub()));
		}
	}
}
