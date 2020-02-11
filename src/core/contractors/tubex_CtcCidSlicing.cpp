/*
 *  CtcBox class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes
 */
#include "tubex_CtcCidSlicing.h"

using namespace std;
using namespace ibex;



namespace tubex
{
	CtcCidSlicing::CtcCidSlicing(int scid, double prec, ibex::Function *fnc): scid(scid), prec(prec), fnc(fnc){
		/*check inputs*/
		assert(scid > 0.);
		assert(prec >= 0);
	}

	void CtcCidSlicing::contract(TubeVector& x, TubeVector& v, TPropagation t_propa, int cid_gate){
		/*check if everything is ok*/
		assert(x.size() == v.size());
		assert(x.domain() == v.domain());
		assert(TubeVector::same_slicing(x, v));
		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;
		for (int i = 0 ; i < x.size() ; i++){
			x_slice.push_back(x[i].first_slice());
			v_slice.push_back(v[i].first_slice());
		}

		/*Defining the sub-contractor Ctc_Derive*/
		/*todo: define these as input?*/
		CtcDeriv ctc_deriv;

		/*for each tube, go all over the slices*/
		while(x_slice[0] != NULL){
			/*iteration step, made for each subslice at each tube x*/
			bool fix_point_n;
			do{
				fix_point_n=false;
				for (int i = 0 ; i < x.size() ; i++){
					std::vector<ibex::Interval> x_subslices;
					x_subslices.clear();
					/*create the slices*/
					create_slices(*x_slice[i],x_subslices, cid_gate);

					/*For each slice on $t$ compute the corresponding the hull */
					Interval hull_input_x = Interval::EMPTY_SET; Interval hull_input_v = Interval::EMPTY_SET;
					Interval hull_output_x = Interval::EMPTY_SET; Interval hull_output_v = Interval::EMPTY_SET;
					Interval hull_codomain_x = Interval::EMPTY_SET; Interval hull_codomain_v = Interval::EMPTY_SET;

					for (int j = 0 ; j < x_subslices.size() ; j++){

						/*Temporal slices on $x$ and $v$*/
						Slice aux_slice_x(*x_slice[i]);
						Slice aux_slice_v(*v_slice[i]);
//
						if (cid_gate == input_gate){
							aux_slice_x.set_input_gate(x_subslices[j]);
						}
						else if (cid_gate == output_gate){
							aux_slice_x.set_output_gate(x_subslices[j]);
						}
						else if (cid_gate == codomain){
							aux_slice_x.set_envelope(x_subslices[j]);
						}

						/*Fixpoint for each sub-slice at each tube*/
						Interval sx;
						do
						{
							sx = aux_slice_x.codomain();
							ctc_deriv.contract(aux_slice_x, aux_slice_v);
							ctc_fwdbwd_slices(aux_slice_x, aux_slice_v, x_slice, v_slice, i);
						} while(std::abs(sx.diam()-aux_slice_x.codomain().diam())>get_prec());


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
				}
				/*is not necessary for 1-dimensional problems*/
				if (x.size() == 1)
					fix_point_n=false;
			} while(fix_point_n);

			if (t_propa & FORWARD){
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i] = x_slice[i]->next_slice();
					v_slice[i] = v_slice[i]->next_slice();
				}
			}
			else if (t_propa & BACKWARD){
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i] = x_slice[i]->prev_slice();
					v_slice[i] = v_slice[i]->prev_slice();
				}
			}
		}
	}


	void CtcCidSlicing::ctc_fwdbwd_slices(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos)
	{

		/*envelope*/
		IntervalVector envelope(x_slice.size());

		for (int i = 0 ; i < x_slice.size() ; i++){
			if (i==pos)
				envelope[i] = x.codomain();
			else
				envelope[i] = x_slice[i]->codomain();
		}
			v.set_envelope(fnc->eval_vector(envelope)[pos]);
	}

	double CtcCidSlicing::get_scid(){
		return this->scid;
	}

	double CtcCidSlicing::get_prec(){
		return this->prec;
	}

	void CtcCidSlicing::create_slices(Slice& x_slice, std::vector<ibex::Interval> & x_slices, int cid_gate){

		/*Varcid in the input gate*/
		if (cid_gate== input_gate){
			double size_interval = x_slice.input_gate().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.input_gate().lb()+i*size_interval,x_slice.input_gate().lb()+size_interval*(i+1)));
			}
		}

		/*Varcid in the output gate*/
		else if (cid_gate == output_gate){
			double size_interval = x_slice.output_gate().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.output_gate().lb()+i*size_interval,x_slice.output_gate().lb()+size_interval*(i+1)));
			}
		}

		/*Varcid in the codomain*/
		else if (cid_gate == codomain){
			double size_interval = x_slice.codomain().diam()/get_scid();
			for (int i = 0 ; i < get_scid() ;i++){
				x_slices.push_back(Interval(x_slice.codomain().lb()+i*size_interval,x_slice.codomain().lb()+size_interval*(i+1)));
			}
		}
	}

	void CtcCidSlicing::change_scid(int scid){
		this->scid = scid;
	}

	void CtcCidSlicing::change_prec(double prec){
		this->prec = prec;
	}
}

