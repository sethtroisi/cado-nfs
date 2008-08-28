#!/usr/bin/perl

# This writes the appropriate Makefile bits that build up a target.

sub usage {
	die "Usage: sub_mkf_maker.pl .<makefile name>.mkf";
}

usage unless scalar @ARGV;

$_ = shift(@ARGV);
my $submkf = $_;
usage unless /^\.(.*)\.mkf$/;

my $target=$1;
my $model;
my $modext="";

if ($target =~ /^(.*?)-([a-z]+)-model$/) {
	$target = $1;
	$model = $2;
	$modext = "-$model";
}

my $sstem;
my $lstem;
my $modext_caps;

($modext_caps = $modext) =~ tr/a-z-/A-Z_/;
($sstem = $target) =~ tr/a-z-/A-Z_/;
$lstem = "$sstem$modext_caps";


open FILE, ">$submkf";
open STDOUT, ">&FILE";

if ($target =~ /^lib/) {
	my $libtarget;
	if (defined($model) && $model eq "so") {
		$libtail=".so";
	} else {
		$libtail="$modext.a";
	}
	print <<EOF;
ifndef	${sstem}_BASE
${sstem}_BASE:=$target
endif
${sstem}_LIBTARGET:=\$(${sstem}_BASE)$libtail
${lstem}_LIB:=\$(BINARY_DIR)\$(${sstem}_LIBTARGET)
LIB_TARGETS+=\$(${lstem}_LIB)
${sstem}_CSOURCES:=\$(filter %.c,\$(${sstem}_SOURCES))
${sstem}_CPPSOURCES:=\$(filter %.cpp,\$(${sstem}_SOURCES))
${sstem}_PCHS:=\$(wildcard \$(${sstem}_CPPSOURCES:.cpp=_pch.hpp))
${lstem}_DEPS:=\$(addprefix \$(DEPDIR),\$(addsuffix ${modext}.d,\$(${sstem}_SOURCES)))
${lstem}_OBJS:=\\
	\$(patsubst %.cpp,\$(BINARY_DIR)%${modext}.o,\$(${sstem}_CPPSOURCES)) \\
	\$(patsubst %.c,\$(BINARY_DIR)%${modext}.o,\$(${sstem}_CSOURCES))
ALLDEPS+=\$(${lstem}_DEPS)
EOF
	if ($model ne 'so') {
		print <<EOF;
\$(${lstem}_LIB):: \$(${lstem}_OBJS) \$(${lstem}_XDEPENDS)
	-rm -f \$@
	\$(call INFORM,LINK,\$(${sstem}_BASE)$libtail\$(PMODEL $model))\$(ARQCV) \$@ \$(${lstem}_OBJS)
EOF
	} else {
		print <<EOF;
\$(${lstem}_LIB):: \$(${lstem}_OBJS) \$(${lstem}_XDEPENDS)
	\$(call INFORM,LINK,\$(${sstem}_BASE)$libtail\$(PMODEL $model))\$(CXX) -shared \$(MODEL_SO_CXXFLAGS) -o \$@ \$(${lstem}_OBJS)
EOF
	}
	print <<EOF;
ifneq (\$(BINARY_DIR),)
.PHONY:	\$(${sstem}_LIBTARGET)
\$(${sstem}_LIBTARGET):
	\@\$(MAKE) --no-print-directory \$(${lstem}_LIB)
endif
CLEAN_PATTERN+=\$(${lstem}_DEPS)
CLEAN_PATTERN+=\$(${lstem}_OBJS) \$(${lstem}_LIB)
${sstem}_GCHS:=\$(${sstem}_PCHS:.hpp=.hpp.gch)
CLEAN_PATTERN+=\$(${sstem}_GCHS) \$(${sstem}_GCHS:%=/tmp/%)
DCLEAN_PATTERN+=$submkf
# vim:ft=make
EOF
} else {
	my $link_addxx = "\$(MODEL__LXXFLAGS)";
	my $link_add = "\$(MODEL__LFLAGS)";
	if ($model) {
		$link_addxx = "\$(MODEL${modext_caps}_LXXFLAGS)";
		$link_add = "\$(MODEL${modext_caps}_LFLAGS)";
	}
	print <<EOF;
ifndef	${sstem}_BASE
${sstem}_BASE:=$target
endif
${sstem}_TARGET:=\$(${sstem}_BASE)$modext
${lstem}_BIN:=\$(BINARY_DIR)\$(${sstem}_TARGET)
BIN_TARGETS+=\$(${lstem}_BIN)
${sstem}_CSOURCES:=\$(filter %.c,\$(${sstem}_SOURCES))
${sstem}_CPPSOURCES:=\$(filter %.cpp,\$(${sstem}_SOURCES))
${sstem}_PCHS:=\$(wildcard \$(${sstem}_CPPSOURCES:.cpp=_pch.hpp))
${lstem}_DEPS:=\$(addprefix \$(DEPDIR),\$(addsuffix ${modext}.d,\$(${sstem}_SOURCES)))
${lstem}_OBJS:=\\
	\$(patsubst %.cpp,\$(BINARY_DIR)%${modext}.o,\$(${sstem}_CPPSOURCES)) \\
	\$(patsubst %.c,\$(BINARY_DIR)%${modext}.o,\$(${sstem}_CSOURCES))
ALLDEPS+=\$(${lstem}_DEPS)
\$(${lstem}_BIN):: \$(${lstem}_OBJS) \$(${lstem}_XDEPENDS)
	\$(call INFORM,LINK,\$(${sstem}_BASE)$modext\$(PMODEL $model))\$(CXX) $link_addxx \$(PFLAGS) \$(OFLAGS) -o \$@ \$(${lstem}_OBJS) \$(LFLAGS) \$(${lstem}_LFLAGS) \$(${lstem}_LIBS)
ifneq (\$(BINARY_DIR),)
.PHONY:	\$(${sstem}_TARGET)
\$(${sstem}_TARGET):
	\@\$(MAKE) --no-print-directory \$(${lstem}_BIN)
endif
CLEAN_PATTERN+=\$(${lstem}_DEPS)
CLEAN_PATTERN+=\$(${lstem}_OBJS) \$(${lstem}_BIN)
${sstem}_GCHS:=\$(${sstem}_PCHS:.hpp=.hpp.gch)
CLEAN_PATTERN+=\$(${sstem}_GCHS) \$(${sstem}_GCHS:%=/tmp/%)
DCLEAN_PATTERN+=$submkf
# vim:ft=make
EOF
}

if (!defined($model)) {
	print <<EOF;
ifneq (\$(${lstem}_EXTRA_MODELS),)
-include \$(patsubst %, .${target}-%-model.mkf, \$(${lstem}_EXTRA_MODELS))
endif
EOF
}

# Care to skip sufficiently many lines so that we do not hit vim's
# modeline part twice.
# vim:ft=perl
