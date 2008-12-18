#include <stddef.h>     /* ptrdiff_t */
#include <stdint.h>
#include "rowset_heap.h"

int heap_index_compare(const struct bucket * a, const struct bucket * b)
{
    return (int32_t) a->i - (int32_t) b->i;
}

int heap_compare(const struct bucket * a, const struct bucket * b)
{
    /* Returns whether b is a better candidate to be on top of the heap
     * IOW: Whether b is lighter than a
     * */
    if (b->room == 0) return 0;
    if (a->room == 0) return 1;
    return b->s < a->s;
}


/* FIXME: libstdc++ is GPL, not LGPL. So we're probably not allow to copy
 * code from there (IANAL). It would be safer to rewrite that.
 */

/* The code below is the glibc heap implementation, as copied from
 * /usr/include/c++/4.1.2/bits/stl_heap.h */
/* {{{ */
void
push_heap0(struct bucket * __first,
        ptrdiff_t hole, ptrdiff_t __topIndex,
        const struct bucket * __value)
{
    ptrdiff_t __parent = (hole - 1) / 2;
    while (hole > __topIndex && heap_compare(__first + __parent, __value))
    {
        *(__first + hole) = *(__first + __parent);
        hole = __parent;
        __parent = (hole - 1) / 2;
    }
    *(__first + hole) = *__value;
}

void
push_heap(struct bucket * __first, struct bucket * __last)
{
    struct bucket ref = __last[-1];
    push_heap0(__first, (__last - __first) - 1, 0, & ref);
}

void adjust_heap(struct bucket * __first, ptrdiff_t hole,
        ptrdiff_t __len, const struct bucket * __value)
{
    const ptrdiff_t __topIndex = hole;
    ptrdiff_t right = 2 * hole + 2;
    while (right < __len)
    {
        if (heap_compare(__first + right,
                    __first + (right - 1)))
            right--;
        *(__first + hole) = *(__first + right);
        hole = right;
        right = 2 * (right + 1);
    }
    if (right == __len)
    {
        *(__first + hole) = *(__first + (right - 1));
        hole = right - 1;
    }
    push_heap0(__first, hole, __topIndex, __value);
}

void pop_heap(struct bucket * __first, struct bucket * __last)
{
    struct bucket ref = __last[-1];
    __last[-1] = __first[0];
    adjust_heap(__first, 0, __last - __first, & ref);
}

void make_heap(struct bucket * __first, struct bucket * __last)
{
    if (__last - __first < 2)
        return;

    const ptrdiff_t __len = __last - __first;
    ptrdiff_t __parent = (__len - 2) / 2;
    for(;;) {
        struct bucket ref = __first[__parent];
        adjust_heap(__first, __parent, __len, &ref);
        if (__parent == 0)
            return;
        __parent--;
    }
}
/* end of heap code */
/* }}} */
