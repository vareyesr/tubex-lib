/** 
 *  CtcPicard class
 * ----------------------------------------------------------------------------
 *  \date       2018
 *  \author     Simon Rohou
 *  \copyright  Copyright 2020 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "tubex_CtcPicard.h"

using namespace std;
using namespace ibex;

#define EPSILON std::numeric_limits<float>::epsilon()

namespace tubex
{
  CtcPicard::CtcPicard(float delta)
    : Ctc(), m_delta(delta)
  {
    assert(delta > 0.);
  }
  
  void CtcPicard::contract(const tubex::Fnc& f, Tube& x, TPropagation t_propa)
  {
    assert(f.nb_vars() == f.image_dim());
    assert(f.nb_vars() == 1 && "scalar case");
    // todo: faster implementation in the scalar case?
    TubeVector x_vect(1, x);
    contract(f, x_vect, t_propa);
    x = x_vect[0];
  }


  void CtcPicard::set_picard_subslices (int nsubslices)
  {
    m_picard_subslices=nsubslices;
  }
  void CtcPicard::contract_picard_tubeslice(const tubex::Fnc& f, TubeVector& x, int & k, TPropagation t_propa )
  {
 /*
          if(!x(k).is_unbounded())
            continue;
	  */
          contract_kth_slices(f, x, k, t_propa);
          // NB: all tube components share the same slicing
          // If the slice stays unbounded after the contraction step,
          // then it is sampled and contracted again.
          if(x(k).is_unbounded() && x[0].slice_domain(k).diam() > x.domain().diam() / m_picard_subslices)
          {

            x.sample(x[0].slice_domain(k).mid()); // all the components of the tube are sampled,
            // and sampling time is selected according to the first slice of one of the components,
            // for instance the first one x[0]

            if (t_propa & FORWARD) k--; // the first subslice will be computed
	    else k+=2;  // the second subslice will be computed
          }
  }

  void CtcPicard::contract_picard_slice(const tubex::Fnc& f, TubeVector& x, int  k, TPropagation t_propa )
  {
    assert (t_propa == FORWARD || t_propa ==BACKWARD);
    Interval initdomain=x[0].slice_domain(k) ;
    int kfinished;
    if (t_propa == FORWARD) kfinished=k+1;
    else kfinished=k-1;
    while (k!= kfinished){
          contract_kth_slices(f, x, k, t_propa);
          // NB: all tube components share the same slicing
          // If the slice stays unbounded after the contraction step,
          // then it is sampled and contracted again.
          if(x(k).is_unbounded())
	     {
	      if (x[0].slice_domain(k).diam() > initdomain.diam() / m_picard_subslices)
		{

		  x.sample(x[0].slice_domain(k).mid()); // all the components of the tube are sampled,
            // and sampling time is selected according to the first slice of one of the components,
            // for instance the first one x[0]

		  if (t_propa == FORWARD) {kfinished++;} // the first subslice will be computed; the final index incremented
		  else {k+=1;}  // the second subslice will be computed
		}
	      else break;  // fail to bound the current slice : stop the algorithm
	     }
	  else if  (t_propa == FORWARD) k++;
	  else k--;
    }
  }


  void CtcPicard::contract(const tubex::Fnc& f, TubeVector& x, TPropagation t_propa)
  {
    assert(f.nb_vars() == f.image_dim());
    assert(f.nb_vars() == x.size());

    if(x.is_empty())
      return;

    if((t_propa & FORWARD) && (t_propa & BACKWARD))
    {
      // todo: select best way according to initial conditions
      contract(f, x, FORWARD);
      contract(f, x, BACKWARD);
    }

    else
    {
      TubeVector *first_slicing;
      if(m_preserve_slicing)
        first_slicing = new TubeVector(x);


      if(t_propa & FORWARD)
      {

        for(int k = 0 ; k < x.nb_slices() ; k++)
        {
	  contract_picard_tubeslice(f,x,k, FORWARD);
        }
      }

      if(t_propa & BACKWARD)
      {
        for(int k = x.nb_slices() - 1 ; k >= 0 ; k--)
        {

	  contract_picard_tubeslice(f,x,k, BACKWARD);

        }
      }

      if(m_preserve_slicing)
      {
        first_slicing->set_empty();
        *first_slicing |= x;
        x = *first_slicing;
        delete first_slicing;
      }
    }
  }

  int CtcPicard::picard_iterations() const
  {
    return m_picard_iterations;
  }

  void CtcPicard::contract_kth_slices(const tubex::Fnc& f,
                                      TubeVector& tube,
                                      int k,
                                      TPropagation t_propa)
  {
    assert(!((t_propa & FORWARD) && (t_propa & BACKWARD)) && "forward/backward case not implemented yet");
    assert(f.nb_vars() == f.image_dim());
    assert(f.nb_vars() == tube.size());
    assert(k >= 0 && k < tube.nb_slices());
    if(tube.is_empty())
      return;

    guess_kth_slices_envelope(f, tube, k, t_propa);
    IntervalVector f_eval = f.eval_vector(k, tube); // computed only once
    if(t_propa & FORWARD)
      for(int i = 0 ; i < tube.size() ; i++)
      {
        Slice *s = tube[i].slice(k);
        s->set_output_gate(s->output_gate()
          & (s->input_gate() + s->domain().diam() * f_eval[i]));
      }

    else if(t_propa & BACKWARD)
      for(int i = 0 ; i < tube.size() ; i++)
      {
        Slice *s = tube[i].slice(k);
        s->set_input_gate(s->input_gate()
          & (s->output_gate() - s->domain().diam() * f_eval[i]));
      }
  }

  void CtcPicard::guess_kth_slices_envelope(const tubex::Fnc& f,
                                            TubeVector& tube,
                                            int k,
                                            TPropagation t_propa)
  {
    assert(!((t_propa & FORWARD) && (t_propa & BACKWARD)) && "forward/backward case not implemented yet");
    assert(f.nb_vars() == f.image_dim());
    assert(f.nb_vars() == tube.size());
    assert(k >= 0 && k < tube.nb_slices());

    if(tube.is_empty())
      return;
    
    float delta = m_delta;
    Interval h, t = tube[0].slice_domain(k);
    IntervalVector initial_x = tube(k), x0(tube.size()), xf(x0);

    if(t_propa & FORWARD)
    {
      x0 = tube(t.lb());
      xf = tube(t.ub());
      h = Interval(0., t.diam());
    }

    else if(t_propa & BACKWARD)
    {
      x0 = tube(t.ub());
      xf = tube(t.lb());
      h = Interval(-t.diam(), 0.);
    }

    IntervalVector x_guess(tube.size()), x_enclosure = x0;
    m_picard_iterations = 0;

    do
    {
      m_picard_iterations++;
      x_guess = x_enclosure;

      for(int i = 0 ; i < x_guess.size() ; i++)
        x_guess[i] = x_guess[i].mid()
                   + delta * (x_guess[i] - x_guess[i].mid())
                   + Interval(-EPSILON,EPSILON); // in case of a degenerate box

      if(f.is_intertemporal())
      {
        // Update needed for further computations
        // that may be related to this slice k
        for(int i = 0 ; i < tube.size() ; i++)
          tube[i].slice(k)->set_envelope(x_guess[i] & initial_x[i]);
        x_enclosure = x0 + h * f.eval_vector(k, tube);
      }

      else // faster evaluation without tube update
      {
        IntervalVector input_box(tube.size() + 1);
        input_box[0] = t;
        input_box.put(1, x_guess & initial_x); // todo: perform the intersection before?
        x_enclosure = x0 + h * f.eval_vector(input_box);
      }

      if(x_enclosure.is_unbounded() || x_enclosure.is_empty() || x_guess.is_empty())
      {
        if(f.is_intertemporal())
          for(int i = 0 ; i < tube.size() ; i++)
            tube[i].slice(k)->set_envelope(initial_x[i]); // coming back to the initial state
        break;
      }
    } while(!x_enclosure.is_interior_subset(x_guess));

    // Setting tube's values
    if(!(x_enclosure.is_unbounded() || x_enclosure.is_empty() || x_guess.is_empty()))
      for(int i = 0 ; i < tube.size() ; i++)
        tube[i].slice(k)->set_envelope(initial_x[i] & x_enclosure[i]);

    if(f.is_intertemporal())
    {
      // Restoring ending gate, contracted by setting the envelope
      for(int i = 0 ; i < tube.size() ; i++)
      {
        Slice *s = tube[i].slice(k);
        if(t_propa & FORWARD)  s->set_output_gate(xf[i]);
        if(t_propa & BACKWARD) s->set_input_gate(xf[i]);
        // todo: ^ check this ^
      }
    }
  }
}
