#ifndef __CPPMETA_H__
#define __CPPMETA_H__



// Concatenating macro, which first expands its arguments.
#define __CAT(x, y) x##y
#define   CAT(x, y) __CAT(x, y)


// Count the number of arguments (up to 4).
#define __NUM_ARGS(x4, x3, x2, x1, n, ...) n
#define   NUM_ARGS(...) __NUM_ARGS(__VA_ARGS__, 4, 3, 2, 1, 0, )


// Apply macro <f> to a list of arguments.
#define __FOR_ALL_ARGS_0(f)
#define __FOR_ALL_ARGS_1(f, a)      f(a, 1)
#define __FOR_ALL_ARGS_2(f, a, ...) f(a, 2), __FOR_ALL_ARGS_1(f, __VA_ARGS__)
#define __FOR_ALL_ARGS_3(f, a, ...) f(a, 3), __FOR_ALL_ARGS_2(f, __VA_ARGS__)
#define __FOR_ALL_ARGS_4(f, a, ...) f(a, 4), __FOR_ALL_ARGS_3(f, __VA_ARGS__)
#define __FOR_ALL_ARGS2(n, f, ...)  __FOR_ALL_ARGS_##n(f, __VA_ARGS__)
#define __FOR_ALL_ARGS1(n, f, ...)  __FOR_ALL_ARGS2(n, f, __VA_ARGS__)
#define   FOR_ALL_ARGS(f, ...)      __FOR_ALL_ARGS1(NUM_ARGS(__VA_ARGS__), \
                                                    f, __VA_ARGS__)


// Return <then> if __<var>_<cond> is defined to "," (a single comma),
// otherwise return <else>.
#define __IF2(x1, x2, ret, ...)     ret
#define __IF1(...)                  __IF2(__VA_ARGS__)
#define   IF(var, cond, then, else) __IF1(CAT(CAT(__, var), CAT(_, cond)), \
                                          then, else, )

// Return <val> if __<var>_<cond> is defined to ", <val>" (a comma followed
// by <val>), otherwise return <default>.
#define SWITCH(var, cond, default)  __IF1(, CAT(CAT(__, var), CAT(_, cond)), \
                                          default, )


// Useful conditions.
#define ___EMPTY       ,

#define __CMP_16_16_EQ ,
#define __CMP_32_32_EQ ,
#define __CMP_64_64_EQ ,
#define __CMP_16_32_LE ,
#define __CMP_16_64_LE ,
#define __CMP_32_64_LE ,
#define   CMP(x, y) CAT(CAT(CMP_, x), CAT(_, y))

#endif  /* __CPPMETA_H__ */
