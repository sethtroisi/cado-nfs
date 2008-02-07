//
// Copyright (C) 2006, 2007, 2008 INRIA (French National Institute for Research
// in Computer Science and Control)
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
//

/**
 * \file    common_funcs.c
 * \author  Jerome Milan
 * \date    Circa February (March?) 2007
 * \version 1.0
 */

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdbool.h>
#include <string.h>

#include "tifa_config.h"
#include "common_funcs.h"

//------------------------------------------------------------------------------
void print_hello_msg(char* algo_name) {
    //
    // Title string
    //
    char title[81];
    snprintf(title, 81, "Integer factorization via the %.40s algorithm", 
             algo_name);
    //
    // An 80-column wide line separator (+1 character for the '\0' of course)
    //
    char linesep[81] = "----------------------------------------"
                       "----------------------------------------";
    //
    // Create the format string to center the title in the line
    //
    char format[8];
    snprintf(format, 81, "%%%lus\n",
             strlen(title) + (strlen(linesep) - strlen(title))/2);

    printf("\n");
    printf("%s\n", linesep);
    printf(format, title);
    printf("%s\n", linesep);
    printf("%s (%s)\n\n", TIFA_FULLNAME, TIFA_VERSION);
}
//------------------------------------------------------------------------------
void print_bye_msg() {
    printf("Bye-bye!\n");
    printf("\n");
}
//------------------------------------------------------------------------------




