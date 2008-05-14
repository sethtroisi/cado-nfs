#ifndef PRECONDITIONER_HPP_
#define PRECONDITIONER_HPP_

#include <stdint.h>
#include <fstream>
#include <string>
#include <vector>
#include <utility>

struct preconditioner_base {
    typedef std::vector<std::pair<uint32_t, int32_t> > srclist_t;
    typedef std::vector<int32_t> dstlist_t;
    typedef std::pair<srclist_t, dstlist_t> combination_t;
    typedef std::vector<combination_t> comblist_t;

    typedef comblist_t::const_iterator clit_t;
    typedef srclist_t::const_iterator sit_t;
    typedef dstlist_t::const_iterator dit_t;

    comblist_t combinations;

    bool init(uint i0, uint i1, std::string const& fn);
};

template<typename traits> class preconditioner : public preconditioner_base
{
    uint i0,i1;
    public:
    bool init(uint xi0, uint xi1, std::string const& fn) {
        i0 = xi0;
        i1 = xi1;
        return preconditioner_base::init(i0, i1, fn);
    }
    void operator()(typename traits::scalar_t * src) const {
        for(clit_t i = combinations.begin() ; i != combinations.end() ; i++) {
            typename traits::wide_scalar_t tmp;
            traits::zero(tmp);
            unsigned int acc = 0;
            for(sit_t j = i->first.begin() ; j != i->first.end() ; j++) {
                traits::addmul(tmp, src[j->first], j->second);
                if (++acc == traits::max_accumulate) {
                    traits::reduce(tmp, tmp);
                    acc = 1;
                }
            }
            for(dit_t j = i->second.begin() ; j != i->second.end() ; j++) {
                typename traits::wide_scalar_t u;
                traits::copy(&u, &tmp, 1);
                traits::addmul(u, src[*j], 1);
                traits::reduce(src[*j], u);
            }
        }
    }
};

/* vim: set sw=4 et sta: */

#endif  /* PRECONDITIONER_HPP_ */
