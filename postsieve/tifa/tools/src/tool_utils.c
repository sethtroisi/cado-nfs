//
// Copyright (C) 2006, 2007 INRIA (French National Institute for Research in
// Computer Science and Control)
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
 * \file    tool_utils.c
 * \author  Jerome Milan
 * \date    Wed Mar 22 2006
 * \version 1.0
 */

#if HAVE_CONFIG_H
#    include <tifa_config.h>
#endif

#include <stdlib.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>

#include "tool_utils.h"
#include "array.h"

//------------------------------------------------------------------------------
bool is_a_number(const char* const str_n, uint32_t length) {
    
    if (0 == isdigit(str_n[0])) {
        return false;
    }
    for (uint32_t i = 1; i < length; i++) {
        if (str_n[i] == '\0') {
            return true;
        }
        if (0 == isdigit(str_n[i])) {
            return false;
        }
    }
    return true;
}
//------------------------------------------------------------------------------
void chomp(char* const str_n, uint32_t length) {
    //
    // _NOTE_: This function is actually quite different from the
    //         Perl builtin chomp function. It's name should probably
    //         be changed to avoid confusion...
    //
    for (uint32_t i = 0; i < length; i++) {
        if (str_n[i] == '\n') {
            str_n[i] = '\0';
            break;
        }
    }
}
//------------------------------------------------------------------------------

