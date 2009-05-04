#ifndef MATMUL_TOY_HPP_
#define MATMUL_TOY_HPP_

#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include "fmt.hpp"
#include "files.hpp"
#include "must_open.hpp"
#include "preconditioner.hpp"
#include "traits.hpp"
#include "matrix.hpp"

template < typename traits > class simple_matmul_toy {
    typedef unsigned int uint;
  public:
    typedef typename traits::scalar_t scalar_t;
    typedef typename traits::wide_scalar_t wide_scalar_t;
    typename traits::representation::matrix_rowset mat;

    uint nr, nc;

    preconditioner < traits > precond;

    scalar_t *v;
    scalar_t *w;

  private:
    wide_scalar_t *scrap;
  public:
    bool is_zero(const scalar_t * vec) {
        for (uint i = 0; i < nc; i++) {
            if (!traits::is_zero(vec[i])) {
                return false;
            }
        }
        return true;
    }

    int nb_nonzero(const scalar_t * vec) {
        int res = 0;
        for (uint i = 0; i < nc; i++) {
            res += !traits::is_zero(vec[i]);
        }
        return res;
    }

    void multiply() {
        precond(v);
        mat.template mul < traits > (scrap, v);
        traits::reduce(w, scrap, 0, nr);
        traits::zero(w + nr, nc - nr);
    }
    void init_common(matrix_stats & stat, const std::string& fname) {
        std::vector<matrix_slice> slices(1);
        stat.need_slices(&slices);

        std::ifstream mtx;
        must_open(mtx, fname);

        stat(mtx);
        mat.fill(mtx, slices[0]);
        mtx.close();

        nr = stat.nr;
        nc = stat.nc;
        v = new scalar_t[nc];
        w = new scalar_t[nc];
        scrap = new wide_scalar_t[nc];
        traits::zero(v, nc);
        traits::zero(w, nc);

        precond.init(0, nr, files::precond);
    }
  public:
    simple_matmul_toy(matrix_stats & stat, const std::string& fname) {
        init_common(stat, fname);
    }
    simple_matmul_toy(const std::string& fname) {
        matrix_stats stat;
        init_common(stat, fname);
    }
    ~simple_matmul_toy() {
        delete[]v;
        delete[]w;
        delete[]scrap;
    }
};


#endif	/* MATMUL_TOY_HPP_ */
