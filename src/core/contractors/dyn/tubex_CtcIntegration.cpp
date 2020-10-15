/*
 *  CtcIntegration Class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#include "tubex_CtcIntegration.h"
#include "tubex_CtcBasic.h"
#include "tubex_CtcCid.h"
#include "tubex_CtcCidGuess.h"

using namespace std;
using namespace ibex;


namespace tubex
{
	CtcIntegration::CtcIntegration(const TFnc& fnc, DynCtc* slice_ctr): fnc(fnc), slice_ctr(slice_ctr),finaltime(-1)
	{

	}

	void CtcIntegration::contract(TubeVector& x, TubeVector& v, double time_dom, TimePropag t_propa, Ptype p_type)
	{

		/*check if everything is ok*/
		assert(x.size() == v.size());
		assert(x.tdomain() == v.tdomain());
		assert(TubeVector::same_slicing(x, v));

		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;

		/*set where to start with the contraction in the dom=[t0,tf]*/
		/*if the contraction is from t=t0*/
		if (time_dom <= x.tdomain().lb()){
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].first_slice());
				v_slice.push_back(v[i].first_slice());
			}
		}

		/*if the contraction is from t=tf*/
		else if (time_dom >= x.tdomain().ub()){
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].last_slice());
				v_slice.push_back(v[i].last_slice());
			}
		}

		/*else the time is inside the interval [t0,tf]*/
		else{
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].slice(time_dom));
				v_slice.push_back(v[i].slice(time_dom));
			}
			if (t_propa & TimePropag::FORWARD){
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i]=x_slice[i]->next_slice();
					v_slice[i]=v_slice[i]->next_slice();
				}
			}
		}

		/*counter for slices*/
		int nb_slices;
		if (t_propa & TimePropag::FORWARD) nb_slices = 0;
		else if (t_propa & TimePropag::BACKWARD) nb_slices = x.nb_slices()-1;

		CtcPicard ctc_picard;
		ctc_picard.set_picard_subslices(10); //todo: setter for nb of subslices
		ctc_picard.preserve_slicing(0);
		//		ctc_picard.preserve_slicing(1);
		/*todo: how to start from any point inside the tube for picard?*/
		if (m_slice_picard_mode){
			if ((time_dom == x.tdomain().lb()) || (time_dom == x.tdomain().ub()))
				m_slice_picard_mode = true;
			else
				m_slice_picard_mode = false;
		}
		std::vector<Interval> idiff_values;

		/*for each tube, go all over the slices*/
		while(x_slice[0] != NULL){

			/*if something is unbounded return*/
			if (m_slice_picard_mode){
				for (unsigned int i = 0 ; i < v_slice.size() ; i++){
					if (v_slice[i]->codomain().is_unbounded()){
						ctc_picard.contract_picard_slice(fnc,x,nb_slices,t_propa);
						v=fnc.eval_vector(x);
						v_slice.clear();
						for (int i = 0 ; i < x.size() ; i++)
							v_slice.push_back(v[i].slice(nb_slices));
						break;
					}
				}
				for (unsigned int i = 0 ; i < x_slice.size() ; i++){ // preserve slicing
				  x_slice[i]=x[i].slice(nb_slices);
				}

			}

			bool contract_slice = true;
			for (unsigned int i = 0 ; i < x_slice.size() ; i++){
				if (v_slice[i]->codomain().is_unbounded()){
					if (m_incremental_mode)
						return;
					else{
						contract_slice = false;
						break;
					}
				}
			}

			if (contract_slice){

				if(p_type & integrodiff){
					idiff_values.push_back(Interval(0.,0.));
					if (contract_idiff(x_slice,v_slice,x,nb_slices,idiff_values,t_propa)){
					  if (t_propa & TimePropag::FORWARD)
							finaltime = x_slice[0]->tdomain().lb();
						else if (t_propa & TimePropag::BACKWARD)
							finaltime = x_slice[0]->tdomain().ub();
						if (m_incremental_mode)
							return;
					}
				}

				else if(dynamic_cast <CtcCid*> (slice_ctr)){
					CtcCid * cid = dynamic_cast <CtcCid*> (slice_ctr);
					if (!cid->contract(x_slice,v_slice,t_propa)){
					  if (t_propa & TimePropag::FORWARD)
							finaltime = x_slice[0]->tdomain().lb();
						else if (t_propa & TimePropag::BACKWARD)
							finaltime = x_slice[0]->tdomain().ub();
						if (m_incremental_mode)
							return;
					}
				}

				else if(dynamic_cast <CtcCidGuess*> (slice_ctr)){
					CtcCidGuess * cidguess = dynamic_cast <CtcCidGuess*> (slice_ctr);
					if (!cidguess->contract(x_slice,v_slice,t_propa)){
					  if (t_propa & TimePropag::FORWARD)
							finaltime = x_slice[0]->tdomain().lb();
					  else if (t_propa & TimePropag::BACKWARD)
							finaltime = x_slice[0]->tdomain().ub();
						if (m_incremental_mode)
							return;
					}
				}
				else if(dynamic_cast <CtcBasic*> (slice_ctr)){
					CtcBasic * basic = dynamic_cast <CtcBasic*> (slice_ctr);
					if (!basic->contract(x_slice,v_slice,t_propa)){
					  if (t_propa & TimePropag::FORWARD)
							finaltime = x_slice[0]->tdomain().lb();
					  else if (t_propa & TimePropag::BACKWARD)
							finaltime = x_slice[0]->tdomain().ub();
						if (m_incremental_mode)
							return;
					}
				}
				else{
					cout << "ERROR: this sub-contractor is not handled by CtcIntegration" << endl;
					return;
				}
			}
			/*continue with the next slice*/
			if (t_propa & TimePropag::FORWARD){
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i] = x_slice[i]->next_slice();
					v_slice[i] = v_slice[i]->next_slice();
				}
			}
			else if (t_propa & TimePropag::BACKWARD){
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i] = x_slice[i]->prev_slice();
					v_slice[i] = v_slice[i]->prev_slice();
				}
			}
			//for picard_slice
			if (t_propa & TimePropag::FORWARD) nb_slices++;
			else if (t_propa & TimePropag::BACKWARD) nb_slices--;
		}


		if (t_propa & TimePropag::FORWARD)
			finaltime = x.tdomain().ub();
		else if (t_propa & TimePropag::BACKWARD)
			finaltime = x.tdomain().lb();
	}

	void CtcIntegration::contract(TubeVector& x, double time_dom, TimePropag t_propa)
	{
		/*v is computed*/
		TubeVector v=x;
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;

		for (int i = 0 ; i < x.size() ; i++){
			x_slice.push_back(x[i].first_slice());
			v_slice.push_back(v[i].first_slice());
		}
		while(x_slice[0] != NULL){
			IntervalVector envelope(x_slice.size()+1);
			envelope[0] = x_slice[0]->tdomain();
			for (unsigned int j = 0 ; j < x_slice.size() ; j++)
				envelope[j+1] = x_slice[j]->codomain();
			envelope = fnc.eval_vector(envelope);
			for (unsigned int j = 0 ; j < x_slice.size() ; j++)
				v_slice[j]->set_envelope(envelope[j]);

			for (int i = 0 ; i < x.size() ; i++)
				x_slice[i] = x_slice[i]->next_slice();
		}

		contract(x,v,time_dom,t_propa);
	}

	double CtcIntegration::get_finaltime()
	{
		return this->finaltime;
	}

	void CtcIntegration::set_incremental_mode(bool incremental_mode){
		this->m_incremental_mode = incremental_mode;
	}

	void CtcIntegration::set_picard_mode(bool slice_picard_mode){
		this->m_slice_picard_mode = slice_picard_mode;
	}

	bool CtcIntegration::contract_idiff(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice,TubeVector x,int id, std::vector<Interval>& idiff_values,TimePropag t_propa){

		Interval to_try(x_slice[0]->tdomain());
		for (unsigned int i = 1 ; i < x_slice.size(); i++)
			assert(to_try == x_slice[i]->tdomain());

		bool fix_point_n;
		bool first_iteration = true;
		CtcDeriv ctc_deriv;
		ctc_deriv.set_fast_mode(true);
		do{
			fix_point_n = false;

			for (unsigned int i = 0 ; i < x_slice.size() ; i++){
//				TubeVector aux_vector_x = x;
				/*Fixpoint for each sub-slice at each tube*/
				double sx;
				/*without polygons*/

				ctc_deriv.contract(*x_slice[i], *v_slice[i],t_propa);
				if (t_propa & TimePropag::FORWARD){
					do{
						sx=x[i].volume();
						Interval aux_codomain;
						/*todo: how to make this general?*/
						Interval integral_value = -5. * x.integral(x_slice[0]->tdomain().lb(),x_slice[0]->tdomain().ub())[i];
						IntervalVector envelope(x_slice.size());
						for (unsigned int j = 0 ; j < x_slice.size() ; j++)
							envelope[j] = x_slice[j]->codomain();
						aux_codomain = fnc.eval_vector(envelope)[i]+ integral_value;
						if (id==0) idiff_values[id]=integral_value;
						else if (id>0) idiff_values[id]=integral_value+idiff_values[id-1];
						if (id > 0 ) aux_codomain+=idiff_values[id-1];
						v_slice[i]->set_envelope(aux_codomain);
						ctc_deriv.contract(*x_slice[i], *v_slice[i],t_propa);
					}while(x[i].volume() != sx);
				}
				if (x_slice[i]->is_empty()) return false;
			}
			if ((first_iteration) && !(fix_point_n))
				return false;

			first_iteration = false;

			/*For 1-dimensional problems it is not necessary to repeat the process*/
			if (x_slice.size() == 1)
				fix_point_n=false;
		} while(fix_point_n);

		return true;
	}
}
