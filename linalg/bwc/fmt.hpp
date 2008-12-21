#ifndef FORMATTER_HPP_
#define FORMATTER_HPP_

#include <iomanip>
#include <string>
#include <sstream>
#include <ostream>
#include <stdexcept>
#include <algorithm>

/* Provides type-safe printf-like formatting for C++.
 *
 * Use as : cout << fmt("Got % shells in % buckets") % nsh % nb << "\n";
 *
 * %% is a literal %
 * %[...] is the extended form for the formatting specifier, which has
 * the following form:
 * (<letter><argument>)*
 *
 * <letter> <argument> <effect>
 *     f      char     setfill
 *     w      uint     setw
 *     .      uint     setprecision
 *
 * Extensions can come easily, of course.
 */

/* boost::format is much, much better than this, but also much,
 * much bigger.
 */
struct fmt_base {
    struct error : public std::runtime_error {
        inline error(const char s[], const std::string& v)
            : std::runtime_error(v + ": " + s) {}
    };
};

class fmt;
class fmt_match;

class fmt : public fmt_base {
    friend class fmt_match;
    std::string s;
    int pos;
    public:
    fmt(const std::string& c, int p = 0) : s(c), pos(p) {}
    fmt(const char * c, int p = 0) : s(c), pos(p) {}
    template<class T>
        fmt operator%(const T& arg) const {
            std::string::size_type s0, s1, s2;
            std::ostringstream m;
            s0 = s.find('%', (std::string::size_type) pos);
            m << s.substr(0, s0);
            for(; s0 != s.npos && s0 + 1 < s.size() && s[s0 + 1] == '%'; ) {
                m << '%';
                s0 += 2;
                if (s0 == s.npos)
                    throw error("too many arguments", s);
                s1 = s.find('%', s0);
                m << s.substr(s0, s1 - s0);
                s0 = s1;
            }
            if (s0 == s.npos)
                throw error("too many arguments", s);
            s1 = s0 + 1;
            s2 = s0;
            std::string fm;
            if (s1 < s.size() && s[s1] == '[') {
                s2 = s.find(']', s1+1);
                if (s2 == s.npos) {
                    throw error("bad pattern", s);
                }
                fm = s.substr(s1+1, s2-s1-1);
            }

            std::istringstream f(fm);
            for(;!f.eof();) {
                char o;
                int w;
                char fl;
                std::streampos oldpos;
                f >> o;
                if (f.eof()) break;
                switch(o) {
                    case 'f':
                        f >> fl;
                        m << std::setfill(fl);
                        break;
                    case 'w':
                        oldpos = f.tellg();
                        f >> w;
                        if (f.tellg() > oldpos)
                            f.clear();
                        m << std::setw(w);
                        break;
                    case  'h':
                        m << std::hex;
                        break;
                    case  'F':
                        m << std::fixed;
                        break;
                    case  '-':
                        m << std::left;
                        break;
                    case '.':
                        oldpos = f.tellg();
                        f >> w;
                        if (f.tellg() > oldpos)
                            f.clear();
                        m << std::setprecision(w);
                        break;
                    default:
                        throw error("bad format", s);
                }
                if (f.fail()) {
                    throw error("bad format", s);
                }
            }
            m << arg;
            m.flush();
            int newpos = m.str().size();
            m << s.substr(s2 + 1);
            return fmt(m.str(), newpos);
        }
    operator std::string() const {
        // do not forget to transform %%'s in the string tail
        //
        std::string v = s;
        std::string::size_type s0;
        for(s0 = pos ; (s0 = v.find("%", s0)) < v.size() ; ) {
            if (v[s0 + 1] != '%')
                throw error("too few arguments", s);
            v = v.erase(++s0, 1);
        }
        return v;
    }
};

/* This does some sort of scanf. It's really not very advanced, has
 * scores of shortcomings and so on. But it works well enough for trivial
 * cases.
 */
class fmt_match : public fmt_base {
    std::string s;      /* format string */
    std::string t;      /* text string */
    bool ok;
    std::string::size_type spos;
    std::string::size_type tpos;
    public:
    fmt_match(const fmt& s, const std::string& t)
        : s(s.s), t(t)
    { ok = true; spos = tpos = 0; }
    fmt_match(const std::string& s, const std::string& t)
        : s(s), t(t)
    { ok = true; spos = tpos = 0; }
    private:
    // try to match up until the next % sign in the format string.
    void eat() {
        std::string::size_type s0;
        for(;;) {
            s0 = s.find('%', (std::string::size_type) spos);
            s0 = std::min(s0, s.size());
            std::string::size_type len = s0 - spos;
            // see if we do have a match that far.
            if (s.compare(spos, len, t, tpos, len) != 0) {
                ok = false;
                return;
            }
            tpos += s0-spos;
            spos = s0;
            if (spos == s.size()) {
                ok = (tpos == t.size());
                return;
            }
            if (s0 != s.npos && s0 + 1 < s.size() && s[s0+1] == '%') {
                if (t[s0 + tpos - spos] == '%') {
                    spos += 2;
                    tpos ++;
                } else {
                    ok = false;
                    return;
                }
            } else {
                break;
            }
        }
    }
    public:
    operator bool() const { return ok; }
    bool eof() const {
        return spos == s.size() || tpos == t.size();
    }
    template<typename T>
        fmt_match operator%(T& arg) {
            eat();
            if (!ok) return *this;
            /* We have some match information */
            std::istringstream m(std::string(t, tpos));
            std::string fm;
            spos++;
            if (spos < s.size() && s[spos] == '[') {
                std::string::size_type s2;
                s2 = s.find(']', spos);
                if (s2 == s.npos) {
                    throw error("bad pattern", s);
                }
                fm = s.substr(spos, s2-spos);
                spos = s2 + 1;
            }
            std::istringstream f(fm);
            for(;!f.eof();) {
                char o;
                int w;
                char fl;
                std::streampos oldpos;
                f >> o;
                if (f.eof()) break;
                switch(o) {
                    case 'f':
                        f >> fl;
                        // m >> std::setfill(fl);
                        break;
                    case 'w':
                        oldpos = f.tellg();
                        f >> w;
                        if (f.tellg() > oldpos)
                            f.clear();
                        // m >> std::setw(w);
                        break;
                    case  'h':
                        m >> std::hex;
                        break;
                    case  'F':
                        m >> std::fixed;
                        break;
                    case  '-':
                        // m >> std::left;
                        break;
                    case '.':
                        oldpos = f.tellg();
                        f >> w;
                        if (f.tellg() > oldpos)
                            f.clear();
                        // m >> std::setprecision(w);
                        break;
                    default:
                        throw error("bad format", s);
                }
                if (f.fail()) {
                    throw error("bad format", s);
                }
            }
            if (!(m >> arg)) {
                ok = false;
                return *this;
            }
            tpos += m.tellg();
            eat();
            return *this;
        }
};

namespace std {
    inline std::ostream& operator<<(std::ostream& o, const fmt& f)
    {
        return o << (std::string) f;
    }
}

#endif	/* FORMATTER_HPP_ */
