
ACLOCAL_AMFLAGS=-I config --install

lib_LTLIBRARIES=libgf2x.la
libgf2x_la_SOURCES=gf2x.c toom.c fft.c
nobase_include_HEADERS=gf2x.h	\
		gf2x/gf2x-impl.h		\
		gf2x/gf2x-small.h

# All these are generated.
nobase_nodist_include_HEADERS=	\
		gf2x/gf2x-config.h		\
		gf2x/gf2x-tfft.h		\
		gf2x/gf2x-thresholds.h		\
		gf2x/gf2x_mul1.h		\
		gf2x/gf2x_mul2.h		\
		gf2x/gf2x_mul3.h		\
		gf2x/gf2x_mul4.h

libgf2x_la_LDFLAGS=-version-info @gf2x_lib_version@

if ENABLE_NTL_CHECKS
    MAYBE_APPS=apps
endif

SUBDIRS=. tests tuning $(MAYBE_APPS)

DIST_SUBDIRS=tests tuning apps

DISTCLEANFILES=gf2x/gf2x-thresholds.h	\
	       gf2x/gf2x_mul1.h		\
	       gf2x/gf2x_mul2.h		\
	       gf2x/gf2x_mul3.h		\
	       gf2x/gf2x_mul4.h

EXTRA_DIST=
EXTRA_DIST+=config/configfsf.guess
EXTRA_DIST+=config/configfsf.sub
EXTRA_DIST+=version.sh.in
EXTRA_DIST+=BUGS
EXTRA_DIST+=patch-wrt-ntl-5.4.2 inject-ntl-5.4.2
EXTRA_DIST+=hardware
dist-hook:
	-find $(distdir) \( -name .svn -type d \) -o -name "*~" -o -name ".#*" \
                | xargs rm -rf

tune-lowlevel tune-toom tune-fft: all
	cd tuning; $(MAKE) $@
