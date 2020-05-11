/*
 *  CtcIntegroDiff
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#include "tubex_CtcDynBasic.h"
#include "tubex_CtcDynCid.h"
#include "tubex_CtcDynCidGuess.h"
#include "tubex_CtcIntegroD.h"


using namespace std;
using namespace ibex;


namespace tubex
{
	CtcIntegroD::CtcIntegroD(tubex::Fnc& fnc, Ctc* slice_ctr): fnc(fnc), slice_ctr(slice_ctr)
	{

	}

	void CtcIntegroD::contract(TubeVector& x, TubeVector& v,TPropagation t_propa){

		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;

		if (t_propa & FORWARD){
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].first_slice());
				v_slice.push_back(v[i].first_slice());
			}
		}
		else if (t_propa & BACKWARD){
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].last_slice());
				v_slice.push_back(v[i].last_slice());
			}
		}
		int slice_id;
		if (t_propa & FORWARD)
			slice_id=0;
		else if (t_propa & BACKWARD)
			slice_id=x.nb_slices()-1;

		while(x_slice[0] != NULL){
//			//contractors (only this one for the moment..)
//			if(dynamic_cast <CtcIntegroDiff*> (slice_ctr)){
//				CtcIntegroDiff * integrodiff = dynamic_cast <CtcIntegroDiff*> (slice_ctr);
//				if (!integrodiff->contract(x_slice,v_slice,x,slice_id,t_propa)){
//
//				}
//			}
			//todo: add incrementality
			if(dynamic_cast <CtcDynBasic*> (slice_ctr)){
				CtcDynBasic * basic = dynamic_cast <CtcDynBasic*> (slice_ctr);
				if (!basic->contract(x_slice,v_slice,x,slice_id,t_propa)){

				}
			}
			//todo: to change to dyncid
			else if(dynamic_cast <CtcDynBasic*> (slice_ctr)){
				CtcDynBasic * basic = dynamic_cast <CtcDynBasic*> (slice_ctr);
				if (!basic->contract(x_slice,v_slice,t_propa)){

				}
			}
			//todo: to change to dyncidguess
			else if(dynamic_cast <CtcDynBasic*> (slice_ctr)){
				CtcDynBasic * basic = dynamic_cast <CtcDynBasic*> (slice_ctr);
				if (!basic->contract(x_slice,v_slice,t_propa)){

				}
			}
			else{
				cout << "ERROR: this sub-contractor is not handled by CtcIntegroD" << endl;
				return;
			}

			if (t_propa & FORWARD){
				slice_id++;
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i] = x_slice[i]->next_slice();
					v_slice[i] = v_slice[i]->next_slice();
				}
			}
			else if (t_propa & BACKWARD){
				slice_id--;
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i] = x_slice[i]->prev_slice();
					v_slice[i] = v_slice[i]->prev_slice();
				}
			}

		}
	}
}
