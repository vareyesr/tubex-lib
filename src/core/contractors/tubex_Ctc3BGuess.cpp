/*
 *  CtcBox class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes
 */
#include "tubex_Ctc3BGuess.h"

using namespace std;
using namespace ibex;



namespace tubex
{
	Ctc3BGuess::Ctc3BGuess(ibex::Fnc& fnc,int bisections, double prec): fnc(fnc), bisections(bisections), prec(prec){
		assert(bisections > 0.);
		assert(prec >= 0);
	}

	void Ctc3BGuess::contract(TubeVector& x, TubeVector& v, TPropagation t_propa, bool m_report){
		/*check if everything is ok*/
		assert(x.size() == v.size());
		assert(x.domain() == v.domain());
		assert(TubeVector::same_slicing(x, v));


		double old_volume = x.volume();
		/*cpu time measurement*/
		clock_t tStart = clock();
		/*init all the tubes*/
		vector<Slice*> x_slice;
		vector<Slice*> v_slice;
		if (t_propa & FORWARD)
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].first_slice());
				v_slice.push_back(v[i].first_slice());
			}
		else if (t_propa & BACKWARD)
			for (int i = 0 ; i < x.size() ; i++){
				x_slice.push_back(x[i].last_slice());
				v_slice.push_back(v[i].last_slice());
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
					create_slices(*x_slice[i],x_subslices, t_propa);

					/*For each slice on $t$ compute the corresponding the hull */
					Interval hull_input_x = Interval::EMPTY_SET; Interval hull_input_v = Interval::EMPTY_SET;
					Interval hull_output_x = Interval::EMPTY_SET; Interval hull_output_v = Interval::EMPTY_SET;
					Interval hull_codomain_x = Interval::EMPTY_SET; Interval hull_codomain_v = Interval::EMPTY_SET;

					/*oracle loop*/
					for (int j = 0 ; j < x_subslices.size() ; j++){

						/*Temporal slices on $x$ and $v$*/
						Slice aux_slice_x(*x_slice[i]);
						Slice aux_slice_v(*v_slice[i]);
//
						if (t_propa & FORWARD)
							aux_slice_x.set_input_gate(x_subslices[j]);

						else if (t_propa & BACKWARD)
							aux_slice_x.set_output_gate(x_subslices[j]);


						/*Fixpoint for each sub-slice at each tube*/
						Interval sx;
//						ctc_deriv.set_fast_mode(true);
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

					/*3B part*/
					if (x.size() == 1){
						x_slice[i]->set_envelope(hull_codomain_x); v_slice[i]->set_envelope(hull_codomain_v);
						x_slice[i]->set_input_gate(hull_input_x); v_slice[i]->set_input_gate(hull_input_v);
						x_slice[i]->set_output_gate(hull_output_x); v_slice[i]->set_output_gate(hull_output_v);
					}

					else{
						Interval c1,c2; // to store the result

						int nb_3b=x_slice[i]->output_gate().diff(hull_output_x,c1,c2);
						c1 = Interval(hull_output_x.ub(),x_slice[i]->output_gate().ub());
						c2 = Interval(x_slice[i]->output_gate().lb(),hull_output_x.lb());
						/*c1 is upper, c2 is lower*/
						bool upper; bool lower;
						/*for c1*/
						if (c1.diam()>0){
							if (oracle_3B(c1,*x_slice[i],*v_slice[i],x_slice,v_slice,ub,i)){
								/*Replacing the old domains with the new ones*/
								x_slice[i]->set_envelope(Interval(x_slice[i]->codomain().lb(),hull_codomain_x.ub()));
								x_slice[i]->set_input_gate(Interval(x_slice[i]->input_gate().lb(),hull_input_x.ub()));
								x_slice[i]->set_output_gate(Interval(x_slice[i]->output_gate().lb(),hull_output_x.ub()));
							}
							else{
								x_slice[i]->set_input_gate(Interval(hull_input_x.lb(),x_slice[i]->input_gate().ub()));
								ctc_deriv.contract(*x_slice[i], *v_slice[i]);
								ctc_fwdbwd_slices(*x_slice[i], *v_slice[i], x_slice, v_slice, i);
								ctc_deriv.contract(*x_slice[i], *v_slice[i]);
							}
						}
						/*for c2*/
						if (c2.diam()>0){
							if (oracle_3B(c2,*x_slice[i],*v_slice[i],x_slice,v_slice,lb,i)){
								/*Replacing the old domains with the new ones*/
								x_slice[i]->set_envelope(Interval(hull_codomain_x.lb(),x_slice[i]->codomain().ub()));
								x_slice[i]->set_input_gate(Interval(hull_input_x.lb(),x_slice[i]->input_gate().ub()));
								x_slice[i]->set_output_gate(Interval(hull_output_x.lb(),x_slice[i]->output_gate().ub()));
							}
							else{
								x_slice[i]->set_input_gate(Interval(hull_input_x.lb(),x_slice[i]->input_gate().ub()));
								ctc_deriv.contract(*x_slice[i], *v_slice[i]);
								ctc_fwdbwd_slices(*x_slice[i], *v_slice[i], x_slice, v_slice, i);
								ctc_deriv.contract(*x_slice[i], *v_slice[i]);
							}
						}
					}
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
//			cout << "slice: " <<x_slice[0]->domain()<< endl;
		}
		if (m_report)
			report(tStart,x,old_volume);
	}


	void Ctc3BGuess::ctc_fwdbwd_slices(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos)
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

	int Ctc3BGuess::get_bisections(){
		return this->bisections;
	}

	double Ctc3BGuess::get_prec(){
		return this->prec;
	}

	void Ctc3BGuess::create_slices(Slice& x_slice, std::vector<ibex::Interval> & x_slices, TPropagation t_propa){

		/*Varcid in the input gate*/
		if (t_propa & FORWARD){
			if (x_slice.input_gate().diam() == 0)
				x_slices.push_back(Interval(x_slice.input_gate().ub()));
			else{
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

	void Ctc3BGuess::change_bisections(int bisections){
		this->bisections = bisections;
	}

	bool Ctc3BGuess::oracle_3B(Interval c1,Slice& xx_slice,Slice& vv_slice, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice,int bound,int pos){

		CtcDeriv ctc_deriv;

		if (bound == ub){
			/*preprocessing: try with c1 first*/
			Slice aux_slice_x(xx_slice);
			Slice aux_slice_v(vv_slice);
			aux_slice_x.set_output_gate(c1);
			Interval sx;
			do
			{
				sx = aux_slice_x.codomain();
				ctc_deriv.contract(aux_slice_x, aux_slice_v);
				ctc_fwdbwd_slices(aux_slice_x, aux_slice_v, x_slice, v_slice, pos);
			} while(std::abs(sx.diam()-aux_slice_x.codomain().diam())>get_prec());

			if (aux_slice_x.is_empty()){
				Interval c3,c4; // to store the result
				int nb_3b=xx_slice.output_gate().diff(c1,c3,c4);
				xx_slice.set_output_gate(c3);
				return true;
			}

			/*the real method*/
			else{
				bool divide = true;
				int iterations = 0; //iterations to 11?
				while((divide) && (iterations < 25)){
					Interval output = Interval(c1.mid(),c1.ub());
					Slice aux_slice_xx(xx_slice);
					Slice aux_slice_vv(vv_slice);
					aux_slice_xx.set_output_gate(output);
					do
					{
						sx = aux_slice_xx.codomain();
						ctc_deriv.contract(aux_slice_xx, aux_slice_vv);
						ctc_fwdbwd_slices(aux_slice_xx, aux_slice_vv, x_slice, v_slice, pos);
					} while(std::abs(sx.diam()-aux_slice_xx.codomain().diam())>get_prec());

					/*we discard the half, update c1 and we continue*/
					if (aux_slice_xx.is_empty()){
						c1 = Interval (c1.lb(),c1.mid());
						iterations++;
					}
					/*if we cant discard, we apply 8 subslices*/
					else{
//						double diam = c1.diam()/8;
//						for (int i = 0 ; i < 8 ; i++){
//							output = Interval(c1.ub()-diam,c1.ub());
//							Slice aux_slice_xx2(xx_slice);
//							Slice aux_slice_vv2(vv_slice);
//							aux_slice_xx2.set_output_gate(output);
//							do
//							{
//								sx = aux_slice_xx2.codomain();
//								ctc_deriv.contract(aux_slice_xx2, aux_slice_vv2);
//								ctc_fwdbwd_slices(aux_slice_xx2, aux_slice_vv2, x_slice, v_slice, pos);
//							} while(std::abs(sx.diam()-aux_slice_xx2.codomain().diam())>get_prec());
//
//							if (aux_slice_xx2.is_empty())
//								c1=Interval(c1.lb(),c1.ub()-diam);
//							else break;
//						}
						divide = false;
					}
				}

			}

		}
		if (bound == lb){
			Slice aux_slice_x(xx_slice);
			Slice aux_slice_v(vv_slice);
			aux_slice_x.set_output_gate(c1);
			Interval sx;
			do
			{
				sx = aux_slice_x.codomain();
				ctc_deriv.contract(aux_slice_x, aux_slice_v);
				ctc_fwdbwd_slices(aux_slice_x, aux_slice_v, x_slice, v_slice, pos);
			} while(std::abs(sx.diam()-aux_slice_x.codomain().diam())>get_prec());

			if (aux_slice_x.is_empty()){
				Interval c3,c4; // to store the result
				int nb_3b=xx_slice.output_gate().diff(c1,c3,c4);
				xx_slice.set_output_gate(c3);
				return true;
			}

			/*the real method*/
			else{
				bool divide = true;
				int iterations = 0; //iterations to 11?
				while((divide) && (iterations < 25)){
					Interval output = Interval(c1.lb(),c1.mid());
					Slice aux_slice_xx(xx_slice);
					Slice aux_slice_vv(vv_slice);
					aux_slice_xx.set_output_gate(output);
					do
					{
						sx = aux_slice_xx.codomain();
						ctc_deriv.contract(aux_slice_xx, aux_slice_vv);
						ctc_fwdbwd_slices(aux_slice_xx, aux_slice_vv, x_slice, v_slice, pos);
					} while(std::abs(sx.diam()-aux_slice_xx.codomain().diam())>get_prec());

					/*we discard the half, update c1 and we continue*/
					if (aux_slice_xx.is_empty()){
						c1 = Interval(c1.mid(),c1.ub());
						iterations++;
					}
					/*if we cant discard, we apply 8 subslices*/
					else{
//						for (int i = 0 ; i < 4 ; i++){
//							output = Interval(c1.ub()-c1.diam(),c1.ub());
//							Slice aux_slice_xx2(xx_slice);
//							Slice aux_slice_vv2(vv_slice);
//							aux_slice_xx2.set_output_gate(output);
//							do
//							{
//								sx = aux_slice_xx2.codomain();
//								ctc_deriv.contract(aux_slice_xx2, aux_slice_vv2);
//								ctc_fwdbwd_slices(aux_slice_xx2, aux_slice_vv2, x_slice, v_slice, pos);
//							} while(std::abs(sx.diam()-aux_slice_xx2.codomain().diam())>get_prec());
//
//							if (aux_slice_xx2.is_empty())
//								c1=Interval(c1.lb(),c1.ub()-c1.diam());
//							else break;
//						}
						divide = false;
					}
				}

			}
		}
		if (bound == ub){
//			ratio=
			xx_slice.set_output_gate(Interval(xx_slice.output_gate().lb(),c1.ub()));
		}
		else{
			xx_slice.set_output_gate(Interval(c1.lb(),xx_slice.output_gate().ub()));
		}
		return false;
	}

	void Ctc3BGuess::report(clock_t tStart,TubeVector& x,double old_volume){

			cout <<endl<< "----------Results for: " <<	dynamic_cast <ibex::Function&>(fnc)<<"----------"<<endl << endl;
			/*CidSlicing does nothing, */
			if (old_volume == x.volume()){
				cout << "\033[1;31mNo contraction made by 3BGuess!\033[0m\n";
				printf("CPU Time spent by 3BGuess: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
			}
			/*CidSlicing contracts the tube*/
			else{
				double doors_size = 0 ;
				int nb_doors = 0;
				for (int i = 0 ; i < x.size() ; i++){
					Slice* x_slice = x[i].first_slice();
					for (int j = 0 ; j < x[i].nb_slices() ; j++){
						doors_size +=x_slice->output_gate().diam();
						nb_doors++;
						x_slice = x_slice->next_slice();
					}
				}
				cout << "\033[1;31mContraction successful!  -  3BGuess\033[0m\n";
				printf("CPU Time spent by 3BGuess: %.3f (s)\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
				printf("Old Volume: %.7f\n", old_volume);
				printf("New Volume: %.7f\n", x.volume());
				printf("Average size of doors: %f\n\n", (double)doors_size/nb_doors);
			}
		}



	void Ctc3BGuess::change_prec(double prec){
		this->prec = prec;
	}
}
