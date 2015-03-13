#include "cado.h"
#include <cstddef>      /* see https://gcc.gnu.org/gcc-4.9/porting_to.html */
#include <string.h>

#include <math.h>

#include <list>
#include <map>
#include <string>
#include <iostream>
#include <sstream>

#include "las-descent-helpers.h"
#include "utils.h"
#include "relation.h"

using namespace std;

struct helper {
    struct tree_label {
        siever_config sc;
        tree_label() { memset(sc, 0, sizeof(siever_config)); }
        tree_label(siever_config_s const& c) { memcpy(sc, &c, sizeof(siever_config)); }
        void set_label(siever_config_s const& c) { memcpy(sc, &c, sizeof(siever_config)); }
        string operator()() const {
            ostringstream os;
            os << sc->bitsize;
            os << '@' << sc->side;
            // char sn[]="ra";
            // os << sn[sc->side];
            return os.str();
        }
        bool operator<(const tree_label& o) const {
            int d = sc->bitsize - o.sc->bitsize;
            if (d) return d<0;
            d = sc->side - o.sc->side;
            return d<0;
        }
    };
    struct tree {
        tree_label label;
        double spent;
        relation_t * rel;
        list<tree *> children;
        tree() { rel = NULL; }
        ~tree() {
            typedef list<tree *>::iterator it_t;
            for(it_t i = children.begin() ; i != children.end() ; i++)
                delete *i;
            children.clear();
            if (rel) {
                relation_clear(rel);
                free(rel);
            }
        }
    };
    list<tree *> forest;
    list<tree *> current;       /* stack of trees */

    helper() {
    }
    ~helper() {
        typedef list<tree *>::iterator it_t;
        for(it_t i = forest.begin() ; i != forest.end() ; i++)
            delete *i;
        forest.clear();
    }

    void new_node(siever_config_srcptr label, int level) {
        ASSERT_ALWAYS(level == (int) current.size());
        tree * kid = new tree;
        kid->label.set_label(*label);
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

    void found_relation(relation_t * rel) {
        relation_t * nrel = (relation_t*) malloc(sizeof(relation_t));
        current.back()->rel = nrel;
        memset(nrel, 0, sizeof(relation_t));
        relation_copy(nrel, rel);
    }
    int depth() {
        return current.size();
    }
    bool is_successful(tree * t) {
        if (!t->rel)
            return 0;
        typedef list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            if (!is_successful(*i))
                return false;
        }
        return true;
    }
    int tree_depth(tree * t) {
        int d = 0;
        typedef list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            d = max(d, 1 + tree_depth(*i));
        }
        return d;
    }
    int tree_weight(tree * t) {
        int w = 1;
        typedef list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            w += tree_weight(*i);
        }
        return w;
    }
    int display_tree(FILE* o, tree * t, string const& prefix) {
        int res = 1;
        fprintf(o, "%s%s [%1.4f]%s\n",
                prefix.c_str(), t->label().c_str(), t->spent,
                t->rel ? "" : " ###");
        string new_prefix = prefix + "  ";
        typedef list<tree *>::iterator it_t;
        for(it_t i = t->children.begin() ; i != t->children.end() ; i++) {
            if (!display_tree(o, *i, new_prefix))
                res = 0;
        }
        if (t->rel == NULL)
            res = 0;
        return res;
    }
    void display_last_tree(FILE * o) {
        bool xs = is_successful(forest.back());
        display_tree(o, forest.back(), xs ? "# ": "# FAILED ");
    }

    struct collected_stats {
        bool ok;
        double t;
        int d;
        int w;
        collected_stats(bool ok, double t, int d, int w) :
            ok(ok), t(t), d(d), w(w) {}
    };
    void display_all_trees(FILE * o) {
        typedef list<tree *>::iterator it_t;
        int total = 0, good = 0;
        for(it_t i = forest.begin() ; i != forest.end() ; i++, total++) {
            bool xs = is_successful(*i);
            good += display_tree(o, *i, xs ? "# ": "# FAILED ");
        }
        fprintf(o, "# Success %d/%d (%1.2f%%)\n", good, total,
                100.0 * (double) good / total);
        list<collected_stats> foo;
        typedef map<tree_label, list<collected_stats> > stats_t;
        stats_t stats;
        typedef stats_t::iterator sit_t;
        for(it_t i = forest.begin() ; i != forest.end() ; i++, total++) {
            collected_stats w(is_successful(*i),
                    (*i)->spent,
                    tree_depth(*i),
                    tree_weight(*i)
                    );
            sit_t si = stats.find((*i)->label);
            if (si == stats.end()) {
                list<collected_stats> v;
                v.push_back(w);
                stats.insert(make_pair((*i)->label, v));
            } else {
                si->second.push_back(w);
            }
        }
        /* We also use this to provide an updated stats file (whose
         * results are based on actual descents, and thus less
         * speculative) */
        list<string> new_hints;
        for(sit_t si = stats.begin() ; si != stats.end() ; si++) {
            /* Collect and print stats for this sq size */
            fprintf(o, "# Stats for %s\n", si->first().c_str());
            int n = si->second.size();
            int nok = 0;
            double t1 = 0;
            double t2 = 0;
            double tmin = 999999;
            double tmax = 0;
            int d1 = 0;
            int dmin = 0;
            int dmax = 0;
            int w1 = 0;
            int wmin = 0;
            int wmax = 0;
            typedef list<collected_stats>::iterator lit_t;
            for(lit_t i = si->second.begin() ; i != si->second.end() ; i++) {
                nok += i->ok;
                t1 += i->t;
                t2 += i->t * i->t;
                d1 += i->d;
                w1 += i->w;
                tmin = min(tmin, i->t);
                tmax = max(tmax, i->t);
                dmin = min(dmin, i->d);
                dmax = max(dmax, i->d);
                wmin = min(wmin, i->w);
                wmax = max(wmax, i->w);
            }
            fprintf(o, "#   success %.2f (%d/%d)\n", (double) nok / n, nok, n);
            fprintf(o, "#   depth avg %.1f, min %d, max %d\n",
                    (double) d1/n, dmin, dmax);
            fprintf(o, "#   weight avg %.1f, min %d, max %d\n",
                    (double) w1/n, wmin, wmax);
            double tt1 = (double) t1/n;
            double tt2 = (double) t2/n;
            fprintf(o, "#   time avg %.3f, sdev %.3f, min %.3f, max %.3f\n",
                    tt1, sqrt(tt2-tt1*tt1), tmin, tmax);

            ostringstream os;
            os << si->first()
                << " " << tt1
                << " " << (double) nok / n 
                << " I=" << si->first.sc->logI;
            for(int i = 0 ; i < 2 ; i++) {
                os << " " << si->first.sc->sides[i]->lim
                   << "," << si->first.sc->sides[i]->lpb
                   << "," << si->first.sc->sides[i]->mfb
                   << "," << si->first.sc->sides[i]->lambda;
            }
            new_hints.push_back(os.str());
        }
        fprintf(o, "# The following data is an _example_ which can be used to provide a hint file\n");
        for(list<string>::iterator i = new_hints.begin() ; i != new_hints.end() ; i++) {
            fprintf(o, "# %s\n", i->c_str());
        }
    }
};



void * las_descent_helper_alloc()
{
    return (void*) new helper;
}
void las_descent_helper_free(void * h)
{
    if (!h) return;
    delete static_cast<helper*>(h);
}


void las_descent_helper_new_node(void * h, siever_config_srcptr label, int level)
{
    if (!h) return;
    static_cast<helper*>(h)->new_node(label, level);
}

void las_descent_helper_found_relation(void * h, relation_t * rel)
{
    if (!h) return;
    static_cast<helper*>(h)->found_relation(rel);
}

void las_descent_helper_done_node(void * h)
{
    if (!h) return;
    static_cast<helper*>(h)->done_node();
}
int las_descent_helper_current_depth(void * h)
{
    if (!h) return 0;
    return static_cast<helper*>(h)->depth();
}

void las_descent_helper_display_last_tree(void * h, FILE * o)
{
    if (!h) return;
    static_cast<helper*>(h)->display_last_tree(o);
}
void las_descent_helper_display_all_trees(void * h, FILE * o)
{
    if (!h) return;
    static_cast<helper*>(h)->display_all_trees(o);
}
