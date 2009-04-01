Name:           mpich2
Version:        1.1b1
Release:        1%{?dist}
Summary:        A Message Passing Interface (MPI) Implementation
Group:          Development/Libraries
License:        BSD
URL:            http://www.mcs.anl.gov/research/projects/mpich2/
Source0:        http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/%{version}/src/%{name}-%{version}.tar.gz
Requires(post): /sbin/ldconfig, /usr/sbin/alternatives
Requires(preun): /usr/sbin/alternatives

%description
MPICH2 is a high-performance and widely portable implementation of the
Message Passing Interface (MPI) standard (both MPI-1 and MPI-2). The
goals of MPICH2 are: (1) to provide an MPI implementation that
efficiently supports different computation and communication platforms
including commodity clusters (desktop systems, shared-memory systems,
multicore architectures), high-speed networks (10 Gigabit Ethernet,
InfiniBand, Myrinet, Quadrics) and proprietary high-end computing systems
(Blue Gene, Cray, SiCortex) and (2) to enable cutting-edge research in
MPI through an easy-to-extend modular framework for other derived
implementations.

# copied from openmpi
%ifarch i386 ppc
%define mode 32
%endif
%ifarch ia64
%define mode 64
%endif
%ifarch s390
%define mode 31
%endif
%ifarch s390x
%define mode 64
%endif
%ifarch x86_64 ppc64
%define mode 64
%endif

%define mpiprefix       /opt/%{name}-%{version}-%{release}.%{_arch}
%define priority        8

%prep
%setup -q

%configure \
	--prefix=%{mpiprefix}%{_prefix}        \
	--libdir=%{mpiprefix}%{_libdir}        \
	--bindir=%{mpiprefix}%{_bindir}        \
	--sbindir=%{mpiprefix}%{_sbindir}        \
	--includedir=%{mpiprefix}%{_includedir}        \
	--sysconfdir=%{mpiprefix}%{_sysconfdir}        \
	--datadir=%{mpiprefix}%{_datadir}        \
	--mandir=%{mpiprefix}%{_mandir}        \
	--with-pm=hydra \
        --enable-sharedlibs=gcc \
        --enable-cxx    \
	--enable-threads=multiple

%build
make

%install
rm -rf ${RPM_BUILD_ROOT}
make install DESTDIR=${RPM_BUILD_ROOT}
ln -s mpicxx ${RPM_BUILD_ROOT}%{mpiprefix}%{_bindir}/mpic++
sed -e '/BUILDROOT/d' -i ${RPM_BUILD_ROOT}%{mpiprefix}%{_sbindir}/mpeuninstall

echo %{mpiprefix}%{_libdir} > ${RPM_BUILD_ROOT}%{mpiprefix}%{_libdir}/%{name}.ld.conf


%clean
[ ! -z "${RPM_BUILD_ROOT}" ] && rm -rf ${RPM_BUILD_ROOT}

%post
alternatives    \
  --install %{_bindir}/mpirun mpi-run   \
                %{mpiprefix}%{_bindir}/mpirun %{priority} \
  --slave %{_bindir}/mpiexec mpi-exec   \
                %{mpiprefix}%{_bindir}/mpiexec \
  --slave %{_mandir}/man1/mpirun.1.gz mpi-run-man \
		%{mpiprefix}%{_mandir}/man1/mpirun.1.gz \
  --slave %{_mandir}/man1/mpiexec.1.gz mpi-exec-man \
		%{mpiprefix}%{_mandir}/man1/mpiexec.1.gz
alternatives    \
  --install %{_sysconfdir}/ld.so.conf.d/mpi%{mode}.conf mpilibs%{mode}  \
                %{mpiprefix}%{_libdir}/%{name}.ld.conf %{priority}
alternatives    \
  --install %{_bindir}/mpicc mpicc \
                %{mpiprefix}%{_bindir}/mpicc %{priority} \
	--slave %{_bindir}/mpic++ mpic++ \
		%{mpiprefix}%{_bindir}/mpic++  \
	--slave %{_bindir}/mpicxx mpicxx \
		%{mpiprefix}%{_bindir}/mpicxx  \
	--slave %{_bindir}/mpiCC mpiCC \
		%{mpiprefix}%{_bindir}/mpicxx  \
	--slave %{_bindir}/mpif77 mpif77 \
		%{mpiprefix}%{_bindir}/mpif77  \
	--slave %{_bindir}/mpif90 mpif90 \
		%{mpiprefix}%{_bindir}/mpif90
/sbin/ldconfig

%preun
alternatives --remove mpicc %{mpiprefix}%{_bindir}/mpicc
alternatives --remove mpilibs%{mode} %{mpiprefix}%{_libdir}/%{name}.ld.conf 
alternatives --remove mpi-run %{mpiprefix}%{_bindir}/mpirun

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root,-)
%{mpiprefix}

%changelog
* Tue Mar 31 2009 Emmanuel Thomé <Emmanuel.Thome@normalesup.org> 1.1b1-1
- bump

* Sun Dec 21 2008 Emmanuel Thomé <Emmanuel.Thome@normalesup.org> 1.1a2-1
- initial pkg.
