#include <octave/oct.h>
#include <octave/ov-struct.h>

DEFUN_DLD (structureSum5, args, ,"Structure sum attempts")
{
  octave_value retval; //If several output values use octave_value_list
  int nargin = args.length ();
  //octave_stdout << "The program has started ";
  if (nargin == 3)
    {
      octave_scalar_map arg0 = args(0).scalar_map_value ();
      octave_scalar_map arg1 = args(1).scalar_map_value ();
      //octave_stdout << "There is 1 argument ";
      if (! error_state)
        {
          std::string alpha = args(2).string_value ();
          //octave_stdout << "alpha ";
          if (! error_state)
            {
            octave_value tmp = arg0.contents (alpha);
            octave_value tmp1 = arg1.contents (alpha);
            //octave_stdout << "tmp = ";

            if (1)
              {
                //tmp = double(tmp);
                //tmp1 = double(tmp1);

                octave_scalar_map st;
                octave_value tmp2 = tmp*tmp1;

                st.assign ("selected", tmp2);

                retval = octave_value (st);
              }
            else
              error("structdemo: struct does not have a field named alpha\n");

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
