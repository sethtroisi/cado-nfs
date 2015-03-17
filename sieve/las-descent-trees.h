#ifndef LAS_DESCENT_DESCENT_TREES_H_
#define LAS_DESCENT_DESCENT_TREES_H_

#include <string>
#include <sstream>
#include <list>
#include <set>
#include <map>
#include <algorithm>    /* max */

#include "relation.h"
#include "las-forwardtypes.h"
#include "las-types.h"

#ifndef __cplusplus
#error "This is C++-only"
#endif

struct descent_tree {
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
        bool operator<(const tree_label& o) const {
            if (pr_cmp()(pr, o.pr)) return true;
            if (pr_cmp()(o.pr, pr)) return false;
            return side < o.side;
        }
    };
    struct tree {
        tree_label label;
        double spent;
        relation rel;
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
    }
    ~descent_tree() {
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

    void found(relation const& rel) {
        current.back()->rel = rel;
        relation_ab ab = rel;
        visited_t::iterator it = visited.find(current.back()->label);
        if (it == visited.end()) {
            visited_t::mapped_type v;
            v.insert(ab);
            visited.insert(std::make_pair(current.back()->label, v));
        } else {
            it->second.insert(ab);
        }
    }
    bool must_avoid(relation const& rel) const {
        relation_ab ab = rel;
        visited_t::const_iterator it = visited.find(current.back()->label);
        if (it == visited.end())
            return false;
        visited_t::mapped_type const& w(it->second);
        visited_t::mapped_type::const_iterator xt = w.find(ab);
        if (xt == w.end())
            return false;
        return true;
    }
    int depth() {
        return current.size();
    }
    bool is_successful(tree * t) {
        if (!t->rel)
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
