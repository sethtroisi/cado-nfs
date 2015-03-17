#ifndef LAS_DESCENT_DESCENT_TREES_H_
#define LAS_DESCENT_DESCENT_TREES_H_

#include <string>
#include <sstream>
#include <list>
#include <algorithm>    /* max */

#include "relation.h"
#include "las-forwardtypes.h"
#include "las-types.h"

#ifndef __cplusplus
#error "This is C++-only"
#endif

struct descent_tree {
    struct tree_label {
        siever_config sc;
        tree_label() { memset(sc, 0, sizeof(siever_config)); }
        tree_label(siever_config_s const& c) { memcpy(sc, &c, sizeof(siever_config)); }
        void set_label(siever_config_s const& c) { memcpy(sc, &c, sizeof(siever_config)); }
        std::string operator()() const {
            std::ostringstream os;
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
        relation rel;
        std::list<tree *> children;
        tree() { }
        ~tree() {
            typedef std::list<tree *>::iterator it_t;
            for(it_t i = children.begin() ; i != children.end() ; i++)
                delete *i;
            children.clear();
        }
    };
    std::list<tree *> forest;
    std::list<tree *> current;       /* stack of trees */

    descent_tree() {
    }
    ~descent_tree() {
        typedef std::list<tree *>::iterator it_t;
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

    void found(relation const& rel) {
        current.back()->rel = rel;
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
