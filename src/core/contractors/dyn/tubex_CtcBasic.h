/*
*  Author(s):    Victor Reyes, Gilles Trombettoni
*  Created:      Mar 6, 2020
*  Updated:      Sep 11, 2020
*/

#ifndef __TUBEX_CTCBASIC_H__
#define __TUBEX_CTCBASIC_H__

#include "tubex_DynCtc.h"
#include "tubex_Slice.h"
#include "tubex_CtcDeriv.h"
#include <vector>

namespace tubex
{
	/**
	* \class CtcBasic
	* \brief \f$\mathcal{C}_{Basic}\f$ that contracts a tube \f$[x](\cdot)\f$ using the contractors
	* \f$\mathcal{C}_{\frac{d}{dt}}\f$ and \f$\mathcal{C}\f$ at the Slice level.
	*/
	class CtcBasic : public DynCtc{

	public:

		/**
		* \brief Creates a contractor object \f$\mathcal{C}_{Basic}\f$
		*
		* \param fnc a tubex function
		* \param prec the relative precision of the method. By default is 0.
		*/
		CtcBasic(const TFnc& fnc, double prec = 0.);
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
		* \brief Contracts a set of abstract domains
		*
		* This method makes the contractor available in the CN framework.
		*
		* \param v_domains vector of Domain pointers
		*/
		void contract(std::vector<Domain*>& v_domains);
		/**
		 * \brief ctc_fwd manages to make an evaluation of the current Slices in order to contract and update v
		 *
		 * \param x the slice
		 * \param v the corresponding evaluation
		 * \param x_slice vector of slices
		 * \param v_slice vector of evaluations
		 * \param pos the dimension
		 *
		 * \todo remove x and v
		 */
		void ctc_fwd(Slice &x, Slice &v, std::vector<Slice*> x_slice, std::vector<Slice*> v_slice, unsigned int pos);
		/**
		 * \brief This method obtains the current precision of the iterative method.
		 */
		double get_prec();
		/**
		 * \brief This method changes the value of the precision.
		 *
		 * \param prec the precision
		 */
		void set_prec(double prec);
		/**
		 * \brief This method returns the current contraction strategy, true: with fixpoint inside the slice ; false: without.
		 */
		bool get_reasoning_slice();
		/**
		 * \brief Sets the contraction strategy.
		 *
		 * \param reasoning_slice a boolean that sets if a fixpoint is used for the contractors.
		 */
		void set_reasoning_slice(bool reasoning_slice = true);

	private:

		bool m_reasoning_slice = true;
		CtcDeriv ctc_deriv;
		const TFnc& fnc;
		double prec;
	};
}

#endif /* SRC_CORE_CONTRACTORS_TUBEX_CTCBASIC_H_ */
