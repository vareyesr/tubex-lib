/*
*  Author(s):    Victor Reyes, Gilles Trombettoni
*  Created:      Mar 6, 2020
*  Updated:      Sep 11, 2020
*/

#ifndef __TUBEX_CTCCIDGUESS_H__
#define __TUBEX_CTCCIDGUESS_H__

#include "tubex_DynCtc.h"
#include "tubex_Slice.h"
#include "tubex_CtcDeriv.h"
#include <vector>
#include <cmath>

namespace tubex
{

	/**
	* \class CtcCidGuess
	* \brief \f$\mathcal{C}_{CidGuess}\f$ contracts a tube \f$[x](\cdot)\f$ using the contractors
	* \f$\mathcal{C}_{\frac{d}{dt}}\f$ and \f$\mathcal{C}\f$ at the Slice level.
	* It uses only the bounds of the input (resp. output) gate on a given dimension of x
	* in order to contract the output (resp. input) gate.
	*/
	class CtcCidGuess : public DynCtc{

	enum {lb,ub};

	public:
		/**
		* \brief Creates a contractor object \f$\mathcal{C}_{CidGuess}\f$
		*
		* \param fnc a tubex function
		* \param prec the relative precision of the method. By default is $0.05$.
		*/
		CtcCidGuess(const TFnc& fnc, double prec=0.05);
		/**
		 * \brief This method performs a contraction for the TubeVector x.
		 * \pre \f$[\mathbf{x}](\cdot)\f$ and \f$[\mathbf{v}](\cdot)\f$ must share the same dimension,
		 * slicing and tdomain.
		 * \param x_slice vector of slices
		 * \param v_slice vector of evaluations
		 * \param t_propa temporal way of propagation, either forward or backward.
		*/
		bool contract(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa);
		/**
		* \brief ctc_fwd manages to make an evaluation of the current Slices in order to
		* contract and update v
		*
		* todo merge both methods, perhaps add it to CtcIntegration?
		*/
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos);
		/**
		* \brief ctc_fwd manages to make an evaluation of the current Slices in order to
		* contract and update v
		*
		* todo merge both methods, perhaps add it to CtcIntegration?
		*/
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice> x_slice, std::vector<Slice> v_slice, int pos);
		/**
		* \brief For each dimension it creates 2 puntual subintervals (the bounds) to be treated
		*
		* \param x_slice the slice to be treated
		* \param slices vector of subintervals (for storing purposes)
		* \param t_propa temporal way of propagation, either forward or backward.
		*/
		void create_bounds(Slice & x_slice, std::vector<ibex::Interval> & slices, TimePropag t_propa);
		/**
		* \brief This method tries to remove a subinterval using the constraint graph
		*
		* \param to_remove interval to be forced
		* \param bound 0: if it is an upper bound , 1: lower bound of the interval
		* \param pos the dimension
		* \param x_slice the current set of slices
		* \param v_slice the current evaluation of the slices
		* \param tt_propa temporal way of propagation, either forward or backward.
		*/
		void var3Bcheck(ibex::Interval to_remove,int bound, int pos ,std::vector<Slice*> & x_slice,std::vector<Slice*> v_slice,TimePropag t_propa);
		/**
		* \brief It Contract the domains of x by using an AC-3 like propagation algorithm.
		*
		* \param x_slice the slices of x.
		* \param v_slice the slices of v.
		* \param t_propa temporal way of propagation, either forward or backward.
		*
		* todo Use Simon's CN propagation?.
		*/
		void FullPropagationEngine(std::vector<Slice> & x_slice, std::vector<Slice> & v_slice, TimePropag t_propa);
		/**
		* \brief It Contract the domains of x by using iteratively \f$\mathcal{C}_{\frac{d}{dt}}\f$ and \f$\mathcal{C}\f$
		*
		* \param x_slice the slices of x.
		* \param v_slice the slices of v.
		* \param t_propa temporal way of propagation, either forward or backward.
		*/
		void AtomicPropagationEngine(std::vector<Slice> & x_slice, std::vector<Slice> & v_slice, TimePropag t_propa);
		/**
		* \brief It creates \f$2^n\f$ points which are used to guess the input (resp. output) gate.
		*
		* \param x_slice the slice to be treated
		* \param points vector of points (for storing purposes)
		* \param t_propa temporal way of propagation, either forward or backward.
		*/
		void create_corners(std::vector<Slice> x_slices, std::vector< std::vector<double> > & points, TimePropag t_propa);
		/**
		* \brief Creates the set of all possible combination of points (Cartesian product)
		* using the bounds of x.
		*
		* \param bounds all the bounds of the variables in x
		*
		* \return A vector of points (also vectors as can be n dimensional)
		*/
		std::vector<std::vector<double>> cart_product (const std::vector<std::vector<double>>& bounds);
		/**
		*  \brief Getter of the current relative precision.
		*
		*  \return the relative precision
		*/
		double get_prec();
		/**
		* \brief Getter for propagation engine type.
		*
		* \return 0: atomic (simple), 1: complete (stronger, but slower)
		*/
		int get_propagation_engine();
		/**
		*  \brief Getter of the maximum number of iterations for the fixpoint. Default is 50.
		*
		*  \return number of iterations
		*/
		bool get_max_it();
		/**
		* \brief Setter for the relative precision.
		* \param prec the precision.
		*/
		void set_prec(double prec);
		/**
		* \brief Setter for the propagation engine.
		* \param engine can be either 0 (simple) or 1 (full)
		*/
		void set_propagation_engine(int engine);
		/**
		* \brief Setter for the maximum number of iterations in the fixpoint.
		*
		* \param max_it number of iterations
		*/
		void set_max_it(bool max_it);
		/**
		* \brief Setter for the subinterval guess.
		*
		* \param s_strategy can be either 0 (for bounds) or 1 (for corners)
		*/
		void set_s_corn(int s_strategy);
		/**
		* \brief Getter for the subinterval guess.
		*
		* \return 0: bounds, 1: corners
		*/
		int get_s_corn();
		/**
		* \brief Setter of the current policy for discarding the guess in the 3BCheck.
		*
		* \param d_policy can be either 0: no policy, 1: dichotomy
		*/
		void set_dpolicy(int d_policy);
		/**
		* \brief Getter of the current policy for the dichotomy part in the 3BCheck.
		*
		* \return 0: no policy, 1: dichotomy
		*/
		int get_dpolicy();

	private:
		const TFnc& fnc;
		double prec;
		CtcDeriv ctc_deriv;
		int engine = 0;  //by default the propagation engine is atomic (faster)
		int s_strategy = 0 ; //by default the propagation engine is atomic (faster)
		bool max_it = false;
		int d_policy = 0; // 0: nothing , 1: small , 2:big
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTC3BGUESS_H_ */
