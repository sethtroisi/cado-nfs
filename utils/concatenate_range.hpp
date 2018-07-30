#ifndef CONCATENATE_RANGE_HPP_
#define CONCATENATE_RANGE_HPP_

namespace concatenate_range_details {

template<typename C1, typename C2>
struct concatenate_range_impl {
    typedef typename C1::value_type T;
    class const_iterator {
        typedef typename C1::const_iterator it1_t;
        typedef typename C2::const_iterator it2_t;
        it1_t cur1, graft;
        it2_t cur2;
        bool passed;
        public:
        const_iterator() = default;
        const_iterator(it1_t a, it1_t graft, it2_t cur2) : cur1(a), graft(graft), cur2(cur2), passed(false) {}
        const_iterator(it2_t a) : cur2(a), passed(true) {}
        const_iterator operator++() {
            if (passed) {
                ++cur2;
            } else {
                ++cur1;
                if (cur1 == graft) {
                    passed = true;
                }
            }
            return *this;
        }
        T const& operator*() { return passed ? *cur2 : *cur1; }
        bool operator!=(const_iterator const & a) const {
            if (passed) {
                return !a.passed || a.cur2 != cur2;
            } else {
                return a.passed || a.cur1 != cur1;
            }
        }
    };
    private:
    friend class const_iterator;
    public:
    const_iterator _begin;
    const_iterator _end;
    const_iterator begin() const { return _begin; }
    const_iterator end() const { return _end; }
    concatenate_range_impl(C1 const & a, C2 const& b) :
        _begin(a.begin(), a.end(), b.begin()),
        _end(b.end())
    {}
};

}

template<typename C1, typename C2>
concatenate_range_details::concatenate_range_impl<C1, C2>
concatenate_range(C1 const & a, C2 const& b) {
    return concatenate_range_details::concatenate_range_impl<C1, C2>(a, b);
}

/* demo. No copy involved. Some overhead at operator++
#include <vector>
#include <list>
#include <iostream>
int main()
{
    std::list<int> v1 { 1,2,3 };
    std::vector<int> v2 { 4,5,6 };
    for(auto & x : concatenate_range(v1, concatenate_range(v2, v1))) {
        std::cout << x << std::endl;
    }
}
*/

#endif	/* CONCATENATE_RANGE_HPP_ */
