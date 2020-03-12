/*
 *  CtcIntegration Class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#include "tubex_CtcIntegration.h"
#include "tubex_CtcDynBasic.h"
#include "tubex_CtcDynCid.h"
#include "tubex_CtcDynCidGuess.h"

using namespace std;
using namespace ibex;


namespace tubex
{
	CtcIntegration::CtcIntegration(ibex::Fnc& fnc, Ctc* slice_ctr): fnc(fnc), slice_ctr(slice_ctr)
	{

	}

	void CtcIntegration::contract(TubeVector& x, TubeVector& v, double time_dom, TPropagation t_propa, bool m_report)
	{
		/*check if everything is ok*/
		assert(x.size() == v.size());
		assert(x.domain() == v.domain());
		assert(TubeVector::same_slicing(x, v));

		/*cpu time measurement*/
		clock_t tStart = clock();

		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;

		/*set where to start with the contraction in the dom=[t0,tf]*/
		/*if the contraction is from t=t0*/
		if (time_dom <= x.domain().lb()){
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].first_slice());
				v_slice.push_back(v[i].first_slice());
			}
		}

		/*if the contraction is from t=tf*/
		else if (time_dom >= x.domain().ub()){
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
			if (t_propa & FORWARD){
				for (int i = 0 ; i < x.size() ; i++){
					x_slice[i]=x_slice[i]->next_slice();
					v_slice[i]=v_slice[i]->next_slice();
				}
			}
		}

		/*for each tube, go all over the slices*/
		while(x_slice[0] != NULL){

			/*todo: call picard to make it bounded*/
			for (int i = 0 ; i < v_slice.size() ; i++){
				if (v_slice[i]->codomain().is_unbounded())
					return;
			}

			if(dynamic_cast <CtcDynCid*> (slice_ctr)){
				CtcDynCid * cid = dynamic_cast <CtcDynCid*> (slice_ctr);
				if (!cid->contract(x_slice,v_slice,t_propa)){
					if (t_propa & FORWARD)
						finaltime = x_slice[0]->domain().lb();
					else if (t_propa & BACKWARD)
						finaltime = x_slice[0]->domain().ub();
					return;
				}
			}

			else if(dynamic_cast <CtcDynCidGuess*> (slice_ctr)){
				CtcDynCidGuess * cidguess = dynamic_cast <CtcDynCidGuess*> (slice_ctr);
				//todo: implement the contraction
			}
			else if(dynamic_cast <CtcDynBasic*> (slice_ctr)){
				CtcDynBasic * basic = dynamic_cast <CtcDynBasic*> (slice_ctr);
				// todo: implement the contraction
			}
			else{
				cout << "ERROR: this sub-contractor is not handled by CtcIntegration" << endl;
				return;
			}

			/*continue with the next slice*/
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
		if (t_propa & FORWARD)
			finaltime = x.domain().ub();
		else if (t_propa & BACKWARD)
			finaltime = x.domain().lb();
	}

	double CtcIntegration::get_finaltime()
	{
		return this->finaltime;
	}

	void CtcIntegration::report(clock_t tStart,TubeVector& x,double old_volume)
	{
//		cout <<endl<< "----------Results for: " <<	dynamic_cast <ibex::Function&>(fnc)<<"----------"<<endl << endl;
//		/*CidSlicing does nothing, */
//		if (old_volume == x.volume()){
//			cout << "\033[1;31mNo contraction made by 3BGuess!\033[0m\n";
//			printf("CPU Time spent by 3BGuess: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//		}
//		/*CidSlicing contracts the tube*/
//		else{
//			double doors_size = 0 ;
//			int nb_doors = 0;
//			for (int i = 0 ; i < x.size() ; i++){
//				Slice* x_slice = x[i].first_slice();
//				for (int j = 0 ; j < x[i].nb_slices() ; j++){
//					doors_size +=x_slice->output_gate().diam();
//					nb_doors++;
//					x_slice = x_slice->next_slice();
//				}
//			}
//			cout << "\033[1;31mContraction successful!  -  3BGuess\033[0m\n";
//			printf("CPU Time spent by 3BGuess: %.5f (s)\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
//			printf("Old Volume: %.7f\n", old_volume);
//			printf("New Volume: %.7f\n", x.volume());
//			printf("Average size of doors: %f\n\n", (double)doors_size/nb_doors);
}


}
