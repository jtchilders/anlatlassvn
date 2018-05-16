#!/usr/bin/env bash
patch -p0 < ../../../patches/buffered_readin_1.patch
patch -p0 < ../../../patches/mpi_aggregate_2_for-unpatched.patch
patch -p0 < ../../../patches/psi_itmin_1.patch
patch -p0 < ../../../patches/downsize_n_speedup_2.patch
patch -p0 < ../../../patches/init_speedup_1.patch
patch -p0 < ../../../patches/mcatnlo_rssplit_1.patch
patch -p0 < ../../../patches/lorentz_function_cache_1.patch
patch -p0 < ../../../patches/slop_reader_removal_1.patch
patch -p0 < ../../../patches/amegic_init_speedup_1.patch
patch -p0 < ../../../patches/amegic_init_speedup_2_1.patch
