#define USE_MARKOWITZ // says it...!

#define MERGE_LEVEL_MAX 256 // maximum level for a merge; such a large value
                            // is only useful when not using BW

/* controls the strategy for Markowitz merge */
#define M_STRATEGY 3 // 0: finish that mergelevel
                     // 1: change if min weight < mergelevel
                     // 2: jump to minimal possible mergelevel
                     // 3: perform one merge, then check for next min weight

