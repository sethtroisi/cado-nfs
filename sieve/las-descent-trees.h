#ifndef LAS_DESCENT_DESCENT_TREES_H_
#define LAS_DESCENT_DESCENT_TREES_H_

#include <string>
#include <sstream>
#include <list>
#include <set>
#include <map>
#include <pthread.h>
#include <algorithm>    /* max */
#include <cmath>        /* isfinite is c99 and std::isfinite is c++11 ;
                         * it's not totally clear that #include <cmath> +
                         * accessing std::isfinite works.
                         */

#include "relation.h"
#include "las-forwardtypes.h"
#include "las-types.h"

#ifndef __cplusplus
#error "This is C++-only"
#endif

struct descent_tree {
    private:
        /* we have const members which need to lock the mutex */
        mutable pthread_mutex_t tree_lock;
    public:
    static double grace_time_ratio;
    struct tree_label {
        int side;
        relation::pr pr;
        tree_label() { }
        tree_label(int side, relation::pr const& pr ) : side(side), pr(pr) {}
        tree_label(int side, mpz_srcptr p, mpz_srcptr r) :side(side), pr(p, r) {}
        std::string operator()() const {
            std::ostringstream os;
            os << mpz_sizeinbase(pr.p, 2) << '@' << side;
            return os.str();
        }
        std::string fullname() const {
            char * str;
            gmp_asprintf(&str, "%d %Zd %Zd", side, pr.p, pr.r);
            std::string s = str;
            free(str);
            return s;
        }
        bool operator<(const tree_label& o) const {
            if (pr_cmp()(pr, o.pr)) return true;
            if (pr_cmp()(o.pr, pr)) return false;
            return side < o.side;
        }
    };
    /* For descent mode: we compute the expected time to finish given the
     * factor sizes, and deduce a deadline.  Assuming that not all
     * encountered factors are below the factor base bound, if we expect
     * an additional time T to finish the decomposition, we keep looking
     * for a better decomposition for a grace time which is computed as
     * x*T, for some configurable ratio x (one might think of x=0.2 for
     * instance. x is the grace_time_ratio member), which defines a
     * ``deadline'' for next step.  [If all factors happen to be smooth,
     * the deadline is immediate, of course.] If within the grace period,
     * a new relation is found, with an earlier implied deadline, the
     * deadline is updated. We say that the "decision is taken" when the
     * deadline passes, and the las machinery is told to decide that it
     * should proceed with the descent, and stop processing the current
     * special-q.
     */
    struct candidate_relation {
        relation rel;
        std::vector<std::pair<int, relation::pr> > outstanding;
        double time_left;
        double deadline;
        // bool marked_taken;      /* false until we take the decision */
        candidate_relation() : time_left(INFINITY), deadline(INFINITY) {} // , marked_taken(false) {}
        candidate_relation& operator=(candidate_relation const& o) {
            /* nothing very fancy, except that we keep the old deadline.
             * */
            rel = o.rel;
            outstanding = o.outstanding;
            time_left = o.time_left;
            if (o.deadline < deadline) deadline = o.deadline;
            return *this;
        }
        bool operator<(candidate_relation const& b) const
        {
            if (!rel) return false;
            if (!b.rel) return true;
            if (std::isfinite(time_left)) { return time_left < b.time_left; }
            return outstanding.size() < b.outstanding.size();
        }
        operator bool() const { return (bool) rel; }
        bool decision_taken() const { return (bool) rel && seconds() >= deadline; }
        void set_time_left(double t) {
            time_left = t;
            deadline = seconds() + grace_time_ratio * t;
        }
    };
    struct tree {
        tree_label label;
        double spent;
        candidate_relation contender;
        std::list<tree *> children;
        tree(tree_label const& label) : label(label) { }
        ~tree() {
            typedef std::list<tree *>::iterator it_t;
            for(it_t i = children.begin() ; i != children.end() ; i++)
                delete *i;
            children.clear();
        }
    };
    std::list<tree *> forest;
    std::list<tree *> current;       /* stack of trees */

    /* This is an ugly temporary hack */
    typedef std::map<
                tree_label,
                std::set<relation_ab>
            > visited_t;
    visited_t visited;


    descent_tree() {
        pthread_mutex_init(&tree_lock, NULL);
    }
    ~descent_tree() {
        pthread_mutex_destroy(&tree_lock);
        typedef std::list<tree *>::iterator it_t;
        for(it_t i = forest.begin() ; i != forest.end() ; i++)
            delete *i;
        forest.clear();
    }

    void new_node(las_todo_entry const * doing) {
        int level = doing->depth;
        ASSERT_ALWAYS(level == (int) current.size());
        tree * kid = new tree(tree_label(doing->side, doing->p, doing->r));
        kid->spent = -seconds();
        if (current.empty()) {
            forest.push_back(kid);
            current.push_back(kid);
        } else {
            current.back()->children.push_back(kid);
            current.push_back(kid);
        }
    }

    void done_node() {
        current.back()->spent += seconds();
        current.pop_back();
    }

    candidate_relation const& current_best_candidate() const {
        /* outside multithreaded context */
        return current.back()->contender;
    }

    /* return true if the decision to go to the next step of the descent
     * should be taken now, and register this relation has having been
     * taken */
    bool must_take_decision() {
        pthread_mutex_lock(&tree_lock);
        bool res = current.back()->contender.decision_taken();
        pthread_mutex_unlock(&tree_lock);
        return res;
    }

    void take_decision() {
        /* must be called outside multithreaded context */
        // current.back()->contender.marked_taken = true;
        relation_ab ab = current.back()->contender.rel;
        visited_t::iterator it = visited.find(current.back()->label);
        if (it == visited.end()) {
            visited_t::mapped_type v;
            v.insert(ab);
            visited.insert(std::make_pair(current.back()->label, v));
        } else {
            it->second.insert(ab);
        }
    }

    /* this returns true if the decision should be taken now */
    bool new_candidate_relation(candidate_relation& newcomer)
    {
        pthread_mutex_lock(&tree_lock);
        candidate_relation & tenant(current.back()->contender);
        if (newcomer < tenant) {
            if (newcomer.outstanding.empty()) {
                verbose_output_print(0, 1, "# [descent] Yiippee, splitting done\n");
            } else if (std::isfinite(tenant.deadline)) {
                /* This implies that newcomer.deadline is also finite */
                double delta = tenant.time_left-newcomer.time_left;
                verbose_output_print(0, 1, "# [descent] Improved ETA by %.2f\n", delta);
            } else if (tenant) {
                /* This implies that we have fewer outstanding
                 * special-q's */
                verbose_output_print(0, 1, "# [descent] Improved number of children to split from %u to %u\n",
                        (unsigned int) tenant.outstanding.size(),
                        (unsigned int) newcomer.outstanding.size());
            }
            tenant = newcomer;
            if (!tenant.outstanding.empty()) {
                verbose_output_print(0, 1, "# [descent] still searching for %.2f\n", tenant.deadline - seconds());
            }
        }
        bool res = tenant.decision_taken();
        pthread_mutex_unlock(&tree_lock);
        return res;
    }

    bool must_avoid(relation const& rel) const {
        relation_ab ab = rel;
        pthread_mutex_lock(&tree_lock);
        visited_t::const_iterator it = visited.find(current.back()->label);
        bool answer = it != visited.end() && it->second.find(ab) != it->second.end();
        pthread_mutex_unlock(&tree_lock);
        return answer;
    }
    int depth() {
        return current.size();
    }
    bool is_successful(tree * t) {
        if (!t->contender.rel)
            return false;
        typedef std::list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            if (!is_successful(*i))
                return false;
        }
        return true;
    }
    int tree_depth(tree * t) {
        int d = 0;
        typedef std::list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            d = std::max(d, 1 + tree_depth(*i));
        }
        return d;
    }
    int tree_weight(tree * t) {
        int w = 1;
        typedef std::list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            w += tree_weight(*i);
        }
        return w;
    }
    int display_tree(FILE* o, tree * t, std::string const& prefix);
    void display_last_tree(FILE * o);

    void display_all_trees(FILE * o);
};
#endif	/* LAS_DESCENT_DESCENT_TREES_H_ */
