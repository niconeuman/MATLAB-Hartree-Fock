#include <octave/oct.h>
#include <octave/ov-struct.h>

DEFUN_DLD (primitiveFactorsSSSS, args, ,
           "This function performs a calculation of all"
           "quantities needed to calculate the Boys Function of"
           "order 0 for all primitive quartets in a shell quartet."
           "It requires the following inputs:"
           "basis_a,basis_b,basis_c,basis_d,Boys_Table."
           "The first 4 are structures, and the 5th is a matrix.")
{
  octave_value retval; //If several output values use octave_value_list
  int nargin = args.length ();

  if (nargin == 1)
  {
      octave_scalar_map arg0 = args(0).scalar_map_value ();
      //octave_scalar_map arg1 = args(1).scalar_map_value ();
      if (! error_state){
          std::string alpha;
          if (! error_state){
            octave_value g1 = arg0.contents (alpha);
            //octave_value g2 = arg1.contents (alpha);

            //octave_value a = g1.contents (3);
            //octave_value b = g2.contents (3);

            retval = g1;// + g2;
          }
          else
            error("First condition passed, but second failed");
      }
      else
          error("An error has occurred in primitiveFactorsSSSS");
  }
  else
    print_usage();
  return retval;
}
