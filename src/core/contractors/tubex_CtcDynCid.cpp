/*
 *  CtcDynCid Class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#include "tubex_CtcDynCid.h"

using namespace std;
using namespace ibex;


namespace tubex
{
	CtcDynCid::CtcDynCid(ibex::Fnc& fnc,int scid, double prec): fnc(fnc), scid(scid), prec(prec)
	{
		/*check inputs*/
		assert(scid > 0.);
		assert(prec >= 0);
	}

	void CtcDynCid::contract(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TPropagation t_propa)
	{
		//todo:make an assert which checks the domain of each slice in x_slice.
		bool fix_point_n;
		bool first_iteration = true;

		do{
			for (int i = 0 ; i < x_slice.size() ; i++){

				std::vector<ibex::Interval> x_subslices;
				x_subslices.clear();

				/*create the sub-slices*/
				create_subslices(*x_slice[i],x_subslices, t_propa);

				/*For each slice on $t$ compute the corresponding the hull */
				Interval hull_input_x = Interval::EMPTY_SET; Interval hull_input_v = Interval::EMPTY_SET;
				Interval hull_output_x = Interval::EMPTY_SET; Interval hull_output_v = Interval::EMPTY_SET;
				Interval hull_codomain_x = Interval::EMPTY_SET; Interval hull_codomain_v = Interval::EMPTY_SET;

				/*treat each subslice, make the hull and then intersect*/
				for (int j = 0 ; j < x_subslices.size() ; j++){

					/*Temporal slices on $x$ and $v$*/
					Slice aux_slice_x(*x_slice[i]);
					Slice aux_slice_v(*v_slice[i]);

					/*If the tube is unbounded, then the algorithm stops*/
					if (v_slice[i]->codomain().is_unbounded())
						return;

					if (t_propa & FORWARD)
						aux_slice_x.set_input_gate(x_subslices[j]);
					else if (t_propa & BACKWARD)
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

					} while(sx - aux_slice_x.volume() > get_prec());

					/*The union of the current Slice is made.*/
					hull_input_x |= aux_slice_x.input_gate(); hull_input_v |= aux_slice_v.input_gate();
					hull_output_x |= aux_slice_x.output_gate(); hull_output_v |= aux_slice_v.output_gate();
					hull_codomain_x |= aux_slice_x.codomain(); hull_codomain_v |= aux_slice_v.codomain();
				}

				double aux_envelope = x_slice[i]->codomain().diam();

				/*Replacing the old domains with the new ones*/
				x_slice[i]->set_envelope(hull_codomain_x); v_slice[i]->set_envelope(hull_codomain_v);
				x_slice[i]->set_input_gate(hull_input_x); v_slice[i]->set_input_gate(hull_input_v);
				x_slice[i]->set_output_gate(hull_output_x); v_slice[i]->set_output_gate(hull_output_v);

				if (aux_envelope > x_slice[i]->codomain().diam())
					fix_point_n = true;
				if ((first_iteration) && !(fix_point_n))
					return;

				first_iteration= false;
			}
			if (x_slice.size() == 1)
				fix_point_n=false;
		} while(fix_point_n);

	}
//	void CtcDynCid::contract(TubeVector& x, TubeVector& v, TPropagation t_propa)
//	{
//		/*check if everything is ok*/
//		assert(x.size() == v.size());
//		assert(x.domain() == v.domain());
//		assert(TubeVector::same_slicing(x, v));
//
//		double old_volume = x.volume();
//
//		/*init all the tubes*/
//		vector<Slice*> x_slice;
//		vector<Slice*> v_slice;
//
//		if (t_propa & FORWARD){
//			for (int i = 0 ; i < x.size() ; i++){
//				x_slice.push_back(x[i].first_slice());
//				v_slice.push_back(v[i].first_slice());
//			}
//		}
//		else if (t_propa & BACKWARD){
//			for (int i = 0 ; i < x.size() ; i++){
//				x_slice.push_back(x[i].last_slice());
//				v_slice.push_back(v[i].last_slice());
//			}
//		}
//
//		/*for each tube, go all over the slices*/
//		while(x_slice[0] != NULL){
//			/*iteration step, made for each subslice at each tube x*/
//			bool fix_point_n;
//			bool first_iteration = true;
//			do{
//				fix_point_n=false;
//				for (int i = 0 ; i < x.size() ; i++){
//					std::vector<ibex::Interval> x_subslices;
//					x_subslices.clear();
//
//					/*create the sub-slices*/
//					create_subslices(*x_slice[i],x_subslices, t_propa);
//					/*For each slice on $t$ compute the corresponding the hull */
//					Interval hull_input_x = Interval::EMPTY_SET; Interval hull_input_v = Interval::EMPTY_SET;
//					Interval hull_output_x = Interval::EMPTY_SET; Interval hull_output_v = Interval::EMPTY_SET;
//					Interval hull_codomain_x = Interval::EMPTY_SET; Interval hull_codomain_v = Interval::EMPTY_SET;
//
//
//					for (int j = 0 ; j < x_subslices.size() ; j++){
//
//						/*Temporal slices on $x$ and $v$*/
//						Slice aux_slice_x(*x_slice[i]);
//						Slice aux_slice_v(*v_slice[i]);
//
//						/*If the tube is unbounded, then the algorithm stops*/
//						if (x_slice[i]->codomain().is_unbounded())
//							return;
//
//						if (t_propa & FORWARD)
//							aux_slice_x.set_input_gate(x_subslices[j]);
//						else if (t_propa & BACKWARD)
//							aux_slice_x.set_output_gate(x_subslices[j]);
//
//						/*Fixpoint for each sub-slice at each tube*/
//						double sx;
//						/*without polygons*/
//						if(m_fast_mode)
//							ctc_deriv.set_fast_mode(true);
//						do
//						{
//							sx = aux_slice_x.volume();
//							ctc_deriv.contract(aux_slice_x, aux_slice_v,t_propa);
//							ctc_fwd(aux_slice_x, aux_slice_v, x_slice, v_slice, i);
//
//						} while(sx - aux_slice_x.volume() > get_prec());
//
//						/*The union of the current Slice is made.*/
//						hull_input_x |= aux_slice_x.input_gate(); hull_input_v |= aux_slice_v.input_gate();
//						hull_output_x |= aux_slice_x.output_gate(); hull_output_v |= aux_slice_v.output_gate();
//						hull_codomain_x |= aux_slice_x.codomain(); hull_codomain_v |= aux_slice_v.codomain();
//					}
//					double aux_envelope = x_slice[i]->codomain().diam();
//
//					/*Replacing the old domains with the new ones*/
//					x_slice[i]->set_envelope(hull_codomain_x); v_slice[i]->set_envelope(hull_codomain_v);
//					x_slice[i]->set_input_gate(hull_input_x); v_slice[i]->set_input_gate(hull_input_v);
//					x_slice[i]->set_output_gate(hull_output_x); v_slice[i]->set_output_gate(hull_output_v);
//
//					if (aux_envelope > x_slice[i]->codomain().diam())
//						fix_point_n = true;
//					if ((first_iteration) && !(fix_point_n))
//						return;
//
//					first_iteration= false;
//				}
//				/*is not necessary for 1-dimensional problems*/
//				if (x.size() == 1)
//					fix_point_n=false;
//			} while(fix_point_n);
//
//
//			/*continue with the next slice*/
//			if (t_propa & FORWARD){
//				for (int i = 0 ; i < x.size() ; i++){
//					x_slice[i] = x_slice[i]->next_slice();
//					v_slice[i] = v_slice[i]->next_slice();
//				}
//			}
//			else if (t_propa & BACKWARD){
//				for (int i = 0 ; i < x.size() ; i++){
//					x_slice[i] = x_slice[i]->prev_slice();
//					v_slice[i] = v_slice[i]->prev_slice();
//				}
//			}
//		}
//	}


	void CtcDynCid::ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos)
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

	double CtcDynCid::get_scid()
	{
		return this->scid;
	}

	double CtcDynCid::get_prec()
	{
		return this->prec;
	}

	void CtcDynCid::set_scid(int scid)
	{
		this->scid = scid;
	}

	void CtcDynCid::set_prec(double prec)
	{
		this->prec = prec;
	}

	void CtcDynCid::create_subslices(Slice& x_slice, std::vector<ibex::Interval> & x_slices, TPropagation t_propa)
	{
		/*Varcid in the input gate*/
		if (t_propa & FORWARD){
			double size_interval = x_slice.input_gate().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.input_gate().lb()+i*size_interval,x_slice.input_gate().lb()+size_interval*(i+1)));
			}
		}

		/*Varcid in the output gate*/
		else if (t_propa & BACKWARD){
			double size_interval = x_slice.output_gate().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.output_gate().lb()+i*size_interval,x_slice.output_gate().lb()+size_interval*(i+1)));
			}
		}
	}
}

