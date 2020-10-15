/*
*  Author(s):    Victor Reyes, Gilles Trombettoni
*  Created:      Mar 6, 2020
*  Updated:      Sep 11, 2020
*/

#ifndef __TUBEX_CTCCID_H__
#define __TUBEX_CTCCID_H__

#include "tubex_DynCtc.h"
#include "tubex_Slice.h"
#include "tubex_CtcDeriv.h"
#include <vector>


namespace tubex
{
	/*
	* \class CtcCid
	* \brief \f$\mathcal{C}_{Cid}\f$ that contracts a tube \f$[x](\cdot)\f$ using the contractors
	* \f$\mathcal{C}_{\frac{d}{dt}}\f$ and \f$\mathcal{C}\f$ at the Slice level.
	* It divides the input (or output) gate on a given dimension of x by several subintervals
	* of equal width.
	*/
	class CtcCid : public DynCtc{

	public:
		/**
		* \brief Creates a contractor object \f$\mathcal{C}_{Cid}\f$
		*
		* \param fnc a tubex function
		* \param scid number of subintervals which are used to divide the gate. By default is $8$.
		* \param prec the relative precision of the method. By default is $0$.
		*/
		CtcCid(const TFnc& fnc,int scid=8, double prec=0.);
		/**
		 * \brief This method performs a contraction for a set of slices on a given domain.
		 * \pre \f$[\mathbf{x}](\cdot)\f$ and \f$[\mathbf{v}](\cdot)\f$ must share the same dimension,
		 * slicing and tdomain.
		 * \param x_slice vector of slices
		 * \param v_slice vector of evaluations
		 * \param t_propa temporal way of propagation, either forward or backward.
		*/
		bool contract(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa);
		/**
		* \brief Contracts a set of abstract domains
		*
		* This method makes the contractor available in the CN framework.
		*
		* \param v_domains vector of Domain pointers
		*/
		void contract(std::vector<Domain*>& v_domains);
		/*
		* \brief Creates a certain number of subintervals to be treated
		*
		* \param x_slice the slice to be treated
		* \param slices vector of subintervals (for storing purposes)
		* \param t_propa temporal way of propagation, either forward or backward.
		*/
		void create_subintervals(Slice & x_slice, std::vector<ibex::Interval> & slices, TimePropag t_propa);
		/**
		* \brief ctc_fwd makes an evaluation of the current Slices using the evolution function in
		* order to contract and update v
		*
		* todo: merge both methods
		*/
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos);
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice> x_slice, std::vector<Slice> v_slice, int pos);
		/**
		* \brief It returns the current number of subintervals.
		*/
		double get_scid();
		/**
		*  \brief It returns the current relative precision.
		*/
		double get_prec();
		/**
		* \brief It returns the current propagation engine type. 0: atomic (simple),
		*  1: complete (stronger, but slower)
		*/
		int get_propagation_engine();
		/**
		* \brief Changes the current number of subintervals.
		* \param scid number of subintervals. Must be an integer
		*/
		void set_scid(int scid);
		/**
		* \brief Changes the propagation engine.
		* \param engine can be either 0 (simple) or 1 (full)
		*/
		void set_propagation_engine(int engine);
		/*
		* \brief Changes the value of the relative precision.
		* \param prec the precision.
		*/
		void set_prec(double prec);
		/*
		* \brief It Contract the domains of $x$ by using an AC-3 like propagation algorithm.
		*
		* \param x_slice the slices of x.
		* \param v_slice the slices of v.
		* \param t_propa temporal way of propagation, either forward or backward.
		*
		* todo: use Simon's CN propagation.
		*/
		void FullPropagationEngine(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa);
		/*
		* \brief It Contract the domains of $x$ by using iteratively \f$\mathcal{C}_{\frac{d}{dt}}\f$ and \f$\mathcal{C}\f$
		*
		* \param x_slice the slices of x.
		* \param v_slice the slices of v.
		* \param t_propa temporal way of propagation, either forward or backward.
		*/
		void AtomicPropagationEngine(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, TimePropag t_propa);
	private:
		const TFnc& fnc;
		int scid;
		double prec;
		CtcDeriv ctc_deriv;
		int engine = 0;  //by default the propagation engine is atomic (faster)
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTCCID_H_ */
