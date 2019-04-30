#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Build supernova
#

VERSION=$(shell git describe --tags --always --dirty)

LIBPY=lib/python

LIBGO=$(shell pwd)/lib/go
LIBGO_BINS=reads/bucket_fastq_by_bc reads/sort_fastq_by_bc
export GOPATH=$(shell pwd)/lib/go:$(shell pwd)/tenkit/lib/go

LIBASSEMBLY=$(shell pwd)/lib/assembly
LIBASSEMBLY_BINS=CP DF ParseBarcodedFastqs FastFastbCount MakeFasta ReScaffold

LIBJEMALLOC=$(shell pwd)/lib/jemalloc
LIBJEMALLOC_BUILD=$(LIBJEMALLOC)/build
LIBJEMALLOC_CONF=$(LIBJEMALLOC_BUILD)/configure
LIBJEMALLOC_VERSION=3.6.0
LIBJEMALLOC_VERSIONED := $(LIBJEMALLOC)/$(LIBJEMALLOC_VERSION)
# for htslib
LIBHTSLIB=$(shell pwd)/lib/htslib
LIBHTSLIB_BUILD=$(LIBHTSLIB)/build
LIBHTSLIB_CONF=$(LIBHTSLIB_BUILD)/configure
LIBHTSLIB_VERSION=1.9
LIBHTSLIB_VERSIONED := $(LIBHTSLIB)/$(LIBHTSLIB_VERSION)
export DEPEND_FLAGS := -static

export CS=

.PHONY: $(LIBASSEMBLY_BINS) $(LIBGO_BINS) test clean jemalloc htslib

#
# Targets for development builds.
#
all: $(LIBJEMALLOC_VERSIONED) $(LIBHTSLIB_VERSIONED) $(LIBASSEMBLY_BINS) gapcloser

$(LIBJEMALLOC_VERSIONED): $(LIBJEMALLOC_CONF)
	make -C $(LIBJEMALLOC_BUILD)
	make -C $(LIBJEMALLOC_BUILD) install_lib install_bin install_include DESTDIR=..

$(LIBJEMALLOC_CONF):
	cd $(LIBJEMALLOC_BUILD); \
	./autogen.sh ;\
	CPPFLAGS="-Wno-unused -Wno-maybe-uninitialized" ./configure --prefix=/$(LIBJEMALLOC_VERSION) --with-jemalloc-prefix=je_

jemalloc: $(LIBJEMALLOC_VERSIONED)

$(LIBHTSLIB_VERSIONED): $(LIBHTSLIB_CONF)
	make -j16 -C $(LIBHTSLIB_BUILD)
	make install -C $(LIBHTSLIB_BUILD) DESTDIR=..

$(LIBHTSLIB_CONF):
	cd $(LIBHTSLIB_BUILD); \
	autoheader; \
	autoconf ;\
	./configure --prefix=/$(LIBHTSLIB_VERSION) --disable-libcurl --disable-lzma --disable-bz2

htslib: $(LIBHTSLIB_VERSIONED)

lib/bin:
	mkdir lib/bin

#$(LIBGO_BINS): lib/bin
#	go install -ldflags "-X $@.__VERSION__ '$(VERSION)'" $@
#	cp lib/go/bin/* lib/bin

$(LIBASSEMBLY_BINS): lib/bin $(LIBASSEMBLY)
	$(MAKE) -j16 -C $(LIBASSEMBLY) $@ STATIC=yes
	cp $(LIBASSEMBLY)/bin/$@ lib/bin

gapcloser: 
	$(MAKE) -C gap_closer
	cp gap_closer/gc lib/bin

test:
	MROPATH=$(PWD)/mro \
	    PATH=$(PATH):$(PWD)/bin \
	    PYTHONPATH=$(PYTHONPATH):$(PWD)/lib/python \
	    nosetests --with-xunit --with-coverage --cover-erase --cover-html --cover-xml --cover-package=stages,tenkit,kitten mro/stages/* lib/python

clean: 
	rm -rf $(LIBGO)/bin
	rm -rf $(LIBGO)/pkg
	rm -rf $(LIBASSEMBLY)/bin
	rm -rf lib/bin
	rm -rf $(LIBJEMALLOC_VERSIONED)
	rm -rf $(LIBJEMALLOC_CONF)
	rm -rf $(LIBHTSLIB_VERSIONED)
	rm -rf $(LIBHTSLIB_CONF)
	make -C $(LIBASSEMBLY) cleanAll
	make -C $(LIBJEMALLOC_BUILD) clean
	make -C $(LIBHTSLIB_BUILD) clean
	make -C gap_closer clean


