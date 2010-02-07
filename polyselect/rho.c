/* Approximation of Dickman's rho function.

Copyright 2010 Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* return an approximation of rho(x) with relative error <= 7.5e-9
   for x <= 10, and an absolute error <= 1.2e-12 for x > 10. */
double
dickman_rho (double x)
{
  if (x < 0)
    {
      fprintf (stderr, "Error in dickman_rho, x < 0\n");
      exit (EXIT_FAILURE);
    }

  if (x <= 1.0)
    return 1.0; /* no error */

  if (x <= 2.0) /* 1 < x <= 2 */
    {
      x = 2 * x - 3.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f1 = dickman_rho.power_series(1, 64)
	 then using Maple:
	 p:=numapprox[minimax](f1,x=-1..1,10,1/f1).
	 It gives a relative error less than 1.2e-9
	 (numapprox[infnorm](f1/p-1,x=-1..1)) */
      return .59453489206592253066+(-.33333334045356381396+(.55555549124525324224e-1+(-.12345536626584334599e-1+(.30864445307049269724e-2+(-.82383544119364408582e-3+(.22867526556051420719e-3+(-.63554707267886054080e-4+(.18727717457043009514e-4+(-.73223168705152723341e-5+.21206262864513086463e-5*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 3.0) /* 2 < x <= 3 */
    {
      x = 2 * x - 5.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f2 = dickman_rho.power_series(2, 64)
	 then using Maple:
	 p:=numapprox[minimax](f2,x=-1..1,10,1/f2).
	 It gives a relative error less than 8.9e-10
	 (numapprox[infnorm](f2/p-1,x=-1..1)) */
      return .13031956190159661490+(-.11890697945147816983+(.45224027972316855319e-1+(-.97335515428611864000e-2+(.20773419306178919493e-2+(-.45596184871832967815e-3+(.10336375145598596368e-3+(-.23948090479838959902e-4+(.58842414590359757174e-5+(-.17744830412711835918e-5+.42395344408760226490e-6*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 4.0) /* 3 < x <= 4 */
    {
      x = 2 * x - 7.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f3 = dickman_rho.power_series(3, 64)
	 then using Maple:
	 p:=numapprox[minimax](f3,x=-1..1,10,1/f3).
	 It gives a relative error less than 7.7e-10
	 (numapprox[infnorm](f3/p-1,x=-1..1)) */
      return .16229593252987041396e-1+(-.18617080355889960242e-1+(.98231465619823138710e-2+(-.30890609038683816355e-2+(.67860233545741575724e-3+(-.13691972815251836202e-3+(.27142792736320271161e-4+(-.54020882619058856352e-5+(.11184852221470961276e-5+(-.26822513121459348050e-6+.53524419961799178404e-7*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 5.0) /* 4 < x <= 5 */
    {
      x = 2 * x - 9.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f4 = dickman_rho.power_series(4, 64)
	 then using Maple:
	 p:=numapprox[minimax](f4,x=-1..1,10,1/f4).
	 It gives a relative error less than 8.0e-10
	 (numapprox[infnorm](f4/p-1,x=-1..1)) */
      return .13701177421087314589e-2+(-.18032881447340571269e-2+(.11344648599460666353e-2+(-.44785452702335769551e-3+(.12312893763747451118e-3+(-.26025855424244979420e-4+(.49439561514862370128e-5+(-.89908455922138763242e-6+(.16408096470458054605e-6+(-.32859809253342863007e-7+.55954809519625267389e-8*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 6.0) /* 5 < x <= 6 */
    {
      x = 2 * x - 11.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f5 = dickman_rho.power_series(5, 64)
	 then using Maple:
	 p:=numapprox[minimax](f5,x=-1..1,10,1/f5).
	 It gives a relative error less than 1.1e-9
	 (numapprox[infnorm](f5/p-1,x=-1..1)) */
      return .86018611204769414547e-4+(-.12455615869107014099e-3+(.87629281627689306328e-4+(-.39688578340310188782e-4+(.12884593599298603240e-4+(-.31758435723373651519e-5+(.63480088103905159169e-6+(-.11347189099214240765e-6+(.19390866486609035589e-7+(-.34493644515794744033e-8+.52005397653635760845e-9*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 7.0) /* 6 < x <= 7 */
    {
      x = 2 * x - 13.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f6 = dickman_rho.power_series(6, 64)
	 then using Maple:
	 p:=numapprox[minimax](f6,x=-1..1,10,1/f6).
	 It gives a relative error less than 1.6e-9
	 (numapprox[infnorm](f6/p-1,x=-1..1)) */
      return .42503555236283394611e-5+(-.66168162622648216578e-5+(.50451140426428793943e-5+(-.25056278110193975867e-5+(.90780066639285274491e-6+(-.25409423435688445924e-6+(.56994106166489666850e-7+(-.10719440462879689904e-7+(.18244728209885020524e-8+(-.30691539036795813350e-9+.42848520313019466802e-10*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 8.0) /* 7 < x <= 8 */
    {
      x = 2 * x - 15.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f7 = dickman_rho.power_series(7, 64)
	 then using Maple:
	 p:=numapprox[minimax](f7,x=-1..1,10,1/f7).
	 It gives a relative error less than 2.6e-9
	 (numapprox[infnorm](f7/p-1,x=-1..1)) */
      return .17178674964103956208e-6+(-.28335703551339176870e-6+(.23000575057461623263e-6+(-.12233608702949464058e-6+(.47877500071532198696e-7+(-.14657787264925438894e-7+(.36368896326425931590e-8+(-.74970784803203153691e-9+(.13390834297490308260e-9+(-.22532330111335516039e-10+.30448478436493554124e-11*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 9.0) /* 8 < x <= 9 */
    {
      x = 2 * x - 17.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f8 = dickman_rho.power_series(8, 64)
	 then using Maple:
	 p:=numapprox[minimax](f8,x=-1..1,10,1/f8).
	 It gives a relative error less than 4.4e-9
	 (numapprox[infnorm](f8/p-1,x=-1..1)) */
      return .58405695885954012372e-8+(-.10105102932563601939e-7+(.86312378192120507636e-8+(-.48483948355429166070e-8+(.20129738917787451004e-8+(-.65800969013673567377e-9+(.17591641970556292848e-9+(-.39380340623747158286e-10+(.75914903051312238224e-11+(-.13345077002021028171e-11+.18138420114766605833e-12*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  if (x <= 10.0) /* 9 < x <= 10 */
    {
      x = 2 * x - 19.0; /* -1 <= x <= 1 */
      /* the following approximation was computed with Sage:
	 f9 = dickman_rho.power_series(9, 64)
	 then using Maple:
	 p:=numapprox[minimax](f9,x=-1..1,10,1/f9).
	 It gives a relative error less than 7.5e-9
	 (numapprox[infnorm](f9/p-1,x=-1..1)) */
      return .17063527516986966493e-9+(-.30739839907061929947e-9+(.27401311519313991860e-9+(-.16103964872970941164e-9+(.70152211692303181302e-10+(-.24143747383189013980e-10+(.68287507297199575436e-11+(-.16282751120203734207e-11+(.33676191519721199042e-12+(-.63207992345353337190e-13+.88821884863349801119e-14*x)*x)*x)*x)*x)*x)*x)*x)*x)*x;
    }

  /* for x > 10, Dickman's function differs from the approximation
     below by at most 1.2e-12 */
  return 0.277017183772596 * pow (x, -x);
}
