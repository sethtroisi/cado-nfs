/*
  include ell_arith.h that defines point_t

  implemented functions:
    - edwards_add (R:output_flag, P:edwards_ext, Q:edwards_ext, output_flag)
        R <- P+Q
        output_flag can be edwards_proj, edwards_ext or montgomery

    - edwards_sub (R:output_flag, P:edwards_ext, Q:edwards_ext, output_flag)
        R <- P-Q
        output_flag can be edwards_proj, edwards_ext or montgomery

    - edwards_dbl (R:output_flag, P:edwards_proj, output_flag)
        R <- 2*P
        output_flag can be edwards_proj, edwards_ext

    - edwards_tpl (R:output_flag, P:edwards_proj, output_flag)
        R <- 3*P
        output_flag can be edwards_proj, edwards_ext

    - edwards_dbladd (R:output_flag, P:edwards_proj, Q:edwards_ext, output_flag)
        R <- 2*P+Q
        output_flag can be edwards_proj, edwards_ext, montgomery
        implemented as
          edwards_dbl (R, P, edwards_ext)
          edwards_add (R, R, Q, output_flag)

    - edwards_dblsub (R:output_flag, P:edwards_proj, Q:edwards_ext, output_flag)
        R <- 2*P-Q
        output_flag can be edwards_proj, edwards_ext, montgomery
        implemented as
          edwards_dbl (R, P, edwards_ext)
          edwards_sub (R, R, Q, output_flag)

    - edwards_tpladd (R:output_flag, P:edwards_proj, Q:edwards_ext, output_flag)
        R <- 3*P+Q
        output_flag can be edwards_proj, edwards_ext, montgomery
        implemented as
          edwards_tpl (R, P, edwards_ext)
          edwards_add (R, R, Q, output_flag)

    - edwards_tplsub (R:output_flag, P:edwards_proj, Q:edwards_ext, output_flag)
        R <- 3*P-Q
        output_flag can be edwards_proj, edwards_ext, montgomery
        implemented as
          edwards_tpl (R, P, edwards_ext)
          edwards_sub (R, R, Q, output_flag)
*/
