/*
 *  CtcIntegration class
 * ----------------------------------------------------------------------------
 * 	\date       2020
 *  \authors  	Victor Reyes, Gilles Trombettoni
 */

#ifndef __TUBEX_CTCINTEGRATION_H__
#define __TUBEX_CTCINTEGRATION_H__

#include "tubex_DynCtc.h"
#include "tubex_Slice.h"
#include "tubex_CtcPicard.h"
#include <vector>
#include <utility>


namespace tubex
{

	enum Ptype
	{
	  ode = 0x01,
	  integrodiff = 0x02,
	  dae = 0x03
	};
	/*
	* \class CtcIntegration
	* \brief \f$\mathcal{C}_{Integration}\f$ manages the contraction of a tube \f$[x](\cdot)\f$ using either
	* \f$\mathcal{C}_{Basic}\f$, \f$\mathcal{C}_{Cid}\f$ or \f$\mathcal{C}_{CidGuess}\f$.
	* Currently, it handles Ordinary Differential Equations (IVPs and BVPs , Integro-Differential
	* Equations and Differential Algebraic Equations.
	*/
	class CtcIntegration : public DynCtc{

	public:
		/**
		* \brief Creates a contractor object \f$\mathcal{C}_{Integration}\f$
		*
		* \param fnc a tubex function
		* \param slice_ctr slice-level contractor: either Basic, Cid or CidGuess.
		*/
		CtcIntegration(const TFnc& fnc, DynCtc* slice_ctr);
		/**
		 * \brief This method performs a contraction for the TubeVector x. (v can be omitted)
		 * \pre \f$[\mathbf{x}](\cdot)\f$ and \f$[\mathbf{v}](\cdot)\f$ must share the same dimension,
		 * slicing and tdomain.
		 * \param x a Tube
		 * \param v a Tube (evaluation of x with the evolution function)
		 * \param time_dom the time where the contraction begins
		 * \param t_propa temporal way of propagation, either forward or backward
		 * \param P_type type of equation, either ode, integrodiff or dae
		*/
		void contract(TubeVector& x, TubeVector& v, double time_dom, TimePropag t_propa, Ptype p_type=ode);
		void contract(TubeVector& x, double time_dom, TimePropag t_propa); //todo: include Ptype here
		/**
		* \brief Contracts a set of abstract domains
		*
		* This method makes the contractor available in the CN framework.
		*
		* \param v_domains vector of Domain pointers
		*/
		void contract(std::vector<Domain*>& v_domains);
		/**
		* \brief ctc_fwd manages to make an evaluation of the current Slices in order to
		* contract and update v
		*
		*/
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, int pos);
		/*
		* \brief This method performs a contraction for integro-diff equations.
		*
		* \param x_slice vector of slices
		* \param v_slices evaluations
		* \param id current id of the slice
		* \param i_diff values of the integral for all the partitions
		* \param t_propa temporal way of propagation, either forward or backward
		*/
		bool contract_idiff(std::vector<Slice*> x_slice, std::vector<Slice*> v_slice,TubeVector x, int id, std::vector<ibex::Interval>& idiff_values,TimePropag t_propa);
		/**
		* \brief returns the final time reached during contraction
		*/
		double get_finaltime();
		/**
		* \brief sets if Picard is used before the contraction (required for infinite tubes)
		*
		* \param slice_picard if True, the Picard theorem is called before contraction in order
		* to bound infinite tubes. By default the number of subdivisions is $8$
		*/
		void set_picard_mode(bool slice_picard = false);
		/**
		* \brief sets if contraction is incremental or not
		*
		* \param incremental_mode if True, then when no contraction is obtained on a given slice the
		* method is finished and the last time $t$ is returned.
		*/
		void set_incremental_mode(bool incremental_mode = true);

	private:
		const TFnc& fnc;
		DynCtc* slice_ctr;
		bool m_incremental_mode = true;
		bool m_slice_picard_mode = false;
		double finaltime;
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTCINTEGRATION_H_ */
