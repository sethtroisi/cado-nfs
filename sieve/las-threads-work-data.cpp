#include "cado.h"
#include "las-types.hpp"
#include "las-threads-work-data.hpp"

void nfs_work::thread_data::side_data::allocate_bucket_region()
{
  /* Allocate memory for each side's bucket region. Our intention is to
   * avoid doing it in case we have no factor base. Not that much because
   * of the spared memory, but rather because it's a useful way to trim
   * the irrelevant parts of the code in that case.
   */
  if (!bucket_region)
  bucket_region = (unsigned char *) contiguous_malloc(BUCKET_REGION + MEMSET_MIN);
}

nfs_work::thread_data::thread_data(thread_data const & o)
    : ws(o.ws), sides(o.sides)
{
    SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);
    memcpy(SS, o.SS, BUCKET_REGION);
}

nfs_work::thread_data::side_data::~side_data()
{
  if (bucket_region) contiguous_free(bucket_region);
  bucket_region = NULL;
}

nfs_work::thread_data::thread_data(nfs_work & ws)
    : ws(ws)
{
  /* Allocate memory for the intermediate sum (only one for both sides) */
  SS = (unsigned char *) contiguous_malloc(BUCKET_REGION);
}

nfs_work::thread_data::~thread_data()
{
    contiguous_free(SS);
}

void nfs_work::thread_data::allocate_bucket_regions()
{
    for (auto & S : sides)
        S.allocate_bucket_region();
}

nfs_work::nfs_work(las_info const & _las)
    : nfs_work(_las, _las.nb_threads + 2)
{}
nfs_work::nfs_work(las_info const & _las, int nr_workspaces)
    : las(_las),
    nr_workspaces(nr_workspaces),
    groups { {nr_workspaces}, {nr_workspaces} },
    th(_las.nb_threads, thread_data(*this))
{ }

nfs_work_cofac::nfs_work_cofac(las_info const& las, sieve_info const & si) :
    las(las),
    sc(si.conf),
    doing(si.doing),
    strategies(si.strategies)
{}

/* Prepare to work on sieving a special-q as described by _si.
   This implies allocating all the memory we need for bucket arrays,
   sieve regions, etc. */
void nfs_work::allocate_buckets(sieve_info const & si, nfs_aux & aux, thread_pool & pool)
{
    /* Always allocate the max number of buckets (i.e., as if we were using the
       max value for J), even if we use a smaller J due to a poor q-lattice
       basis */
    /* Take some margin depending on parameters */
    /* Multithreading perturbates the fill-in ratio */
    bkmult_specifier const & multiplier = * si.bk_multiplier;
    verbose_output_print(0, 2, "# Reserving buckets with a multiplier of %s\n",
            multiplier.print_all().c_str());

    bool do_resieve = si.conf.sides[0].lim && si.conf.sides[1].lim;

    for (unsigned int side = 0; side < 2; side++) {
        sieve_info::side_info const & sis(si.sides[side]);
        if (sis.fb->empty()) continue;
        groups[side].allocate_buckets(si.nb_buckets,
                multiplier,
                si.sides[side].fbs->stats.weight,
                si.conf.logI,
                aux, pool, do_resieve);
    }
    pool.drain_queue(2);
}

void nfs_work::allocate_bucket_regions() {
    for(auto & T : th)
        T.allocate_bucket_regions();
}

template <int LEVEL, typename HINT>
double
nfs_work::buckets_max_full()
{
    /* find the most full bucket across all buckets in the bucket array */
    double maxfull_ratio = 0;
    int maxfull_side = -1;
    unsigned int maxfull_index = 0;
    size_t maxfull_updates = 0;
    size_t maxfull_room = 0;
    typedef bucket_array_t<LEVEL, HINT> BA_t;
    for(int side = 0 ; side < 2 ; side++) {
        for (auto const & BA : bucket_arrays<LEVEL, HINT>(side)) {
            unsigned int index;
            const double ratio = BA.max_full(&index);
            if (ratio > maxfull_ratio) {
                maxfull_ratio = ratio;
                maxfull_side = side;
                maxfull_index = index;
                maxfull_updates = BA.nb_of_updates(index);
                maxfull_room = BA.room_allocated_for_updates(index);
            }
        }
    }
    if (maxfull_ratio > 1) {
        int side = maxfull_side;
        auto const & BAs = bucket_arrays<LEVEL, HINT>(side);
        std::ostringstream os;
        os << "bucket " << maxfull_index << " on side " << maxfull_side << ":";
        for (auto const & BA : bucket_arrays<LEVEL, HINT>(side)) {
            os << " " << BA.nb_of_updates(maxfull_index);
        }
        os << " /" << BAs[0].room_allocated_for_updates(maxfull_index);

        auto k = bkmult_specifier::getkey<typename BA_t::update_t>();
        verbose_output_print(0, 1, "# Error: %s buckets are full, worst is %s\n",

                bkmult_specifier::printkey(k).c_str(),
                os.str().c_str());

        throw buckets_are_full(
                k,
                maxfull_index,
                maxfull_updates,
                maxfull_room);
    }
    return maxfull_ratio;
}
double nfs_work::check_buckets_max_full()
{
    double mf0 = 0, mf;
    mf = buckets_max_full<3, shorthint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, shorthint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, longhint_t>();  if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, shorthint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, longhint_t>();  if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<3, emptyhint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, emptyhint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, logphint_t>();  if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, emptyhint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, logphint_t>();  if (mf > mf0) mf0 = mf;
    return mf0;
}

template <typename HINT>
double nfs_work::check_buckets_max_full(int level)
{
    switch(level) {
        case 3:
            if (std::is_same<HINT, shorthint_t>::value)
                return buckets_max_full<3, shorthint_t>();
            else if (std::is_same<HINT, emptyhint_t>::value)
                return buckets_max_full<3, emptyhint_t>();
            else
                ASSERT_ALWAYS(0);
        case 2: return buckets_max_full<2, HINT>();
        case 1: return buckets_max_full<1, HINT>();
    }
    return 0;
}

template double nfs_work::buckets_max_full<1, shorthint_t>();
template double nfs_work::buckets_max_full<2, shorthint_t>();
template double nfs_work::buckets_max_full<3, shorthint_t>();
template double nfs_work::buckets_max_full<1, longhint_t>();
template double nfs_work::buckets_max_full<2, longhint_t>();
template double nfs_work::buckets_max_full<1, emptyhint_t>();
template double nfs_work::buckets_max_full<2, emptyhint_t>();
template double nfs_work::buckets_max_full<3, emptyhint_t>();
template double nfs_work::buckets_max_full<1, logphint_t>();
template double nfs_work::buckets_max_full<2, logphint_t>();
template double nfs_work::check_buckets_max_full<shorthint_t>(int);
template double nfs_work::check_buckets_max_full<longhint_t>(int);
template double nfs_work::check_buckets_max_full<emptyhint_t>(int);
template double nfs_work::check_buckets_max_full<logphint_t>(int);

template <int LEVEL, typename HINT>
void
nfs_work::reset_all_pointers(int side) {
    groups[side].get<LEVEL, HINT>().reset_all_pointers();
}

template void nfs_work::reset_all_pointers<1, shorthint_t>(int);
template void nfs_work::reset_all_pointers<2, shorthint_t>(int);
template void nfs_work::reset_all_pointers<3, shorthint_t>(int);
template void nfs_work::reset_all_pointers<1, emptyhint_t>(int);
template void nfs_work::reset_all_pointers<2, emptyhint_t>(int);
template void nfs_work::reset_all_pointers<3, emptyhint_t>(int);
