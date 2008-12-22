//#define USE_MARKOWITZ // says it...!

#define MERGE_LEVEL_MAX 256 // maximum level for a merge; such a large value
                            // is only useful when not using BW

#define M_STRATEGY 3 // 0: finish that mergelevel
                     // 1: change if min weight < mergelevel
                     // 2: jump to minimal possible mergelevel
                     // 3: perform one merge, then check for next min weight

// 0 means dummy filtering using slow methods
// 1 means:
//  * fast with optimized-though-memory-consuming data structure;
//  * heavy columns are killed.
// 2 means:
//  * fast as above;
//  * heavy columns are present, but not in S[]; this yields more accurate
//    weights.
// 3 means:
//  * reintroduce columns whose weight is <= mergelevelmax (not suggested
//    for large numbers).
#define USE_MERGE_FAST 2

