# Generic part, of the testing Makefiles.
# Possible calls are:
#  make in: generate all input dirs.
#  make failed-check: run check again on all directories with missing report
#                     or failed report.
#  make X.in: generate input dir for directory X.
#  make X.check: generate a report for directory X (if not already existing).
#  make X.recheck: force the creation of the report in directory X.
#  make X.clean: clean the given directory X.
#  make X.diff: make the difference between the output and the reference (with DIFF envvar)
#  make X.updateref: update the reference with the output (prompt the overwrite)

#write here a portable way to run only few test in the Makefile.
#taken as the third solution of the interesting webpage http://gallium.inria.fr/blog/portable-conditionals-in-makefiles/
#thanks to this trick we can test only few test by typing the command
# make check checkonly_that=C O2-Spin etc. CHECK_MODE=custom
#if CHECK_MODE has not been defined promote it to long check
#the following variable is true if we are inside a bzr branch, false otherwise
CHECK_MODE_INTERNAL = $(shell if test -f $(top_srcdir)/branchfile; then echo true; else echo false;fi)

#this is the policy to be given in the case of explicit CHECK_MODE
checkonlyfoo_short_true= short
checkonlyfoo_short_false= short
checkonlyfoo_long_true= long
checkonlyfoo_long_false= long
checkonlyfoo_custom_true= that
checkonlyfoo_custom_false= that
#this is what would happen if the CHECK_MODE variable is not defined
checkonlyfoo__true= long
checkonlyfoo__false= short

checkonly_short=$(SHORT_TESTDIRS)
checkonly_long=$(LONG_TESTDIRS)
#this fixes the default value, if the CHECK_MODE is badly set
checkonly_=$(SHORT_TESTDIRS)


TESTDIRS := ${checkonly_${checkonlyfoo_${CHECK_MODE}_${CHECK_MODE_INTERNAL}}}


# here one might also reset the values for future use, but automake would complain
#checkonlyfoo_short=
#checkonlyfoo_long=
#checkonly_short=
#checkonly_long=
#checkonly_=


#$(TESTDIRS)

ALLDIRS = $(EXTRA_TESTDIRS) $(LONG_TESTDIRS)

INS = $(ALLDIRS:=.in)
RUNS = $(ALLDIRS:=.run)
CHECKS = $(TESTDIRS:=.check) $(TESTDIRS:=.yaml-check)
EXTRA_CHECKS = $(EXTRA_TESTDIRS:=.check) $(EXTRA_TESTDIRS:=.yaml-check)
DIFFS = $(ALLDIRS:=.diff)
UPDATES = $(ALLDIRS:=.updateref)
FAILEDCHECKS = $(TESTDIRS:=.recheck)
CLEANS = $(ALLDIRS:=.clean)

EXTRA_DIST += $(ALLDIRS)

in: $(INS)

check: $(CHECKS) report

complete-check: $(EXTRA_CHECKS) check

diff: $(DIFFS)

update-references: $(UPDATES)

clean: $(CLEANS)

distclean: $(CLEANS)
	rm -rf Makefile

failed-check: $(FAILEDCHECKS) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

#Binary dependencies
%.NEB.out: $(abs_top_builddir)/spred/src/NEB NEB_include.sh NEB_driver.sh
	rm -f neb.it*
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_serial) $(abs_top_builddir)/spred/src/NEB | tee $@
	cat neb.NEB.0*/log.yaml | grep -v "Unable to read mpd.hosts" > log.yaml
	echo "---" >> log.yaml
	grep ":" NEB.NEB.out | grep -v "<BigDFT>" >> log.yaml
	grep -v ":" NEB.NEB.out > tmp-neb.out
	mv tmp-neb.out NEB.NEB.out
	rm -rf neb.NEB.0*
	rm -f gen_output_file velocities_file
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.splsad.out: $(abs_top_builddir)/spred/src/splsad
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/spred/src/splsad -l yes > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.minhop.out: $(abs_top_builddir)/spred/src/global
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/spred/src/global -l yes > $@
#	mv log-mdinput.yaml log.yaml
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"
%.mhgps.out: $(abs_top_builddir)/spred/src/mhgps
	if test -n "${LD_LIBRARY_PATH}" ; then export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; fi ; \
	$(run_parallel) $(abs_top_builddir)/spred/src/mhgps -l yes > $@
	name=`basename $@ .out` ; \
	$(MAKE) -f ../Makefile $$name".post-out"

%.clean:
	@name=`basename $@ .clean` ; dir=$$name-test ; \
	rm -rf $$dir ; rm -f $$name.* ; \
    echo "Clean $$dir"

%.post-in: ;
%.psp: ;
%.post-clean: ;
%.post-out: ;

in_message:
	@if test -n "$(run_ocl)" ; then \
	  echo "==============================================" ; \
	  echo "Will generate a 'input.perf' file to force OCL" ; \
	    if test -n "$(ocl_platform)" ; then \
	      echo "Forcing use of $(ocl_platform)" ; \
	    fi ; \
	    if test -n "$(ocl_devices)" ; then \
	      echo "Forcing use of $(ocl_devices)" ; \
	    fi ; \
	  echo "==============================================" ; \
	fi

$(INS): in_message
	@name=`basename $@ .in` ; dir=$$name-test ; \
	if test ! -d $$dir ; then mkdir $$dir ; fi ; \
	for i in $(srcdir)/$$name/* ; do cp -rf $$i $$dir ; done ; \
	chmod u+w $$dir/* ; \
	if test -n "$(accel_in_message)" -a -n "$(run_ocl)" ; then \
		if test "$(run_ocl)" = "CPU" ; then \
				echo "accel OCLCPU" > $$dir/check.perf ; \
		elif test "$(run_ocl)" = "ACC" ; then \
				echo "accel OCLACC" > $$dir/check.perf ; \
		else \
				echo "accel OCLGPU" > $$dir/check.perf ; \
		fi ; \
		if test -n "$(ocl_platform)" ; then \
				echo "OCL_PLATFORM $(ocl_platform)" >> $$dir/check.perf ; \
		fi ; \
		if test -n "$(ocl_devices)" ; then \
				echo "OCL_DEVICES $(ocl_devices)" >> $$dir/check.perf ; \
		fi ; \
		for i in $$dir/*.out.ref.yaml ; do \
			base=`basename $$i .out.ref.yaml | $(SED) "s/[^_]*_\?\(.*\)$$/\1/"` ; \
	    	 	if test -n "$$base" ; then cat $$dir/check.perf >> $$dir/$$base.perf ; \
			else cat $$dir/check.perf >> $$dir/input.perf ; \
			fi ; \
		done ; \
	fi ; \
	cd $$dir && $(MAKE) -f ../Makefile $$name".psp"; \
	$(MAKE) -f ../Makefile $$dir".post-in"; \
	echo "Input prepared in \"$$dir\" directory, make $$name.run available"
	touch $@


run_message:
	@if test -n "$(run_parallel)" ; then \
	  echo "==============================================" ; \
	  echo "Will run tests in parallel with '$$run_parallel'" ; \
	  echo "==============================================" ; \
	fi

%.run: %.in run_message
	@name=`basename $@ .run` ; dir=$$name-test ; \
        runs="$(srcdir)/$$name/*.ref" ; \
	tgts=`for r in $$runs ; do echo $$(basename $$r .ref)".out"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts ; \
        echo "Tests have run in \"$$dir\" directory, make $$name.check available"
	@touch $@

%.check: %.run %.yaml-check
	@name=`basename $@ .check` ; dir=$$name-test ; \
        chks="$(srcdir)/$$name/*.ref" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref)".report"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@

%.yaml-check: %.run
	@name=`basename $@ .yaml-check` ; dir=$$name-test ; \
        chks="$(srcdir)/$$name/*.ref.yaml" ; \
	tgts=`for c in $$chks ; do echo $$(basename $$c .ref.yaml)".report.yaml"; done` ; \
        cd $$dir && $(MAKE) -f ../Makefile $$tgts
	touch $@


%.diff	: %.run
	@if test -z "$$DIFF" ; then echo "The environment variable DIFF is missing!"; else \
			name=`basename $@ .diff` ; dir=$$name-test ; \
				chks="$(srcdir)/$$name/*.ref" ; \
			for c in $$chks ; do \
			    echo "$$DIFF $$c $$dir/$$(basename $$c .ref).out"; \
				$$DIFF $$c $$dir/$$(basename $$c .ref)".out"; \
			done ; \
				ychks="$(srcdir)/$$name/*.ref.yaml" ; \
			for c in $$ychks ; do base=`basename $$c .out.ref.yaml | $(SED) s/.out// | $(SED) s/.xabs// | $(SED) "s/[^_]*_\?\(.*\)$$/\1/"` ; \
			if test -n "$$base" ; then \
			echo "$$DIFF $$c $$dir/log-$$base.yaml" ; \
			$$DIFF $$c $$dir/log-$$base.yaml; \
			else \
			echo "$$DIFF $$c $$dir/log.yaml" ; \
			$$DIFF $$c $$dir/log.yaml; \
			fi ;\
			done ; \
    fi ; \
	touch $@

%.updateref: #%.run %.diff
	@name=`basename $@ .updateref` ; dir=$$name-test ; \
        chks="$(srcdir)/$$name/*.ref" ; \
	for c in $$chks ; do echo "Update reference with " $$dir/$$(basename $$c .ref)".out"; \
	                     cp -vi $$dir/$$(basename $$c .ref)".out"  $$c;\
	done ; \
        ychks="$(srcdir)/$$name/*.ref.yaml" ; \
	for c in $$ychks ; do base=`basename $$c .out.ref.yaml | $(SED) s/.out// | $(SED) s/.xabs// | $(SED) "s/[^_]*_\?\(.*\)$$/\1/"`  ;\
	if test -n "$$base" ; then \
	echo "Update reference with " $$dir/log-$$base.yaml; \
	                     cp -vi $$dir/log-$$base.yaml $$c;\
	else \
	echo "Update reference with " $$dir/log.yaml; \
	                     cp -vi $$dir/log.yaml $$c;\
	fi ;\
	done ; \
	touch $@

%.recheck: %.in
	@name=`basename $@ .recheck` ; dir=$$name-test ; \
	refs="$$dir/*.ref" ; \
	for r in $$refs ; do \
	  rep=`basename $$r .ref`".report" ; \
	  if ! grep -qs "succeeded\|passed" $$dir/$$rep ; then \
	    target=` basename $$r .ref` ; \
	    rm -f $$dir/$$target".out" $$dir/$$target".report" ; \
	    cd $$dir && $(MAKE) -f ../Makefile $$target".out" $$target".report" && cd - ; \
	  fi \
	done
	touch $*".check"

# Avoid copying in dist the builddir files.
distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | $(SED) 's/[].[^$$\\*]/\\\\&/g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | $(SED) 's/[].[^$$\\*]/\\\\&/g'`; \
	list='$(DISTFILES)'; \
	  dist_files=`for file in $$list; do echo $$file; done | \
	  $(SED) -e "s|^$$srcdirstrip/||;t" \
	      -e "s|^$$topsrcdirstrip/|$(top_builddir)/|;t"`; \
	case $$dist_files in \
	  */*) $(MKDIR_P) `echo "$$dist_files" | \
			   $(SED) '/\//!d;s|^|$(distdir)/|;s,/[^/]*$$,,' | \
			   sort -u` ;; \
	esac; \
	for file in $$dist_files; do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    dir=`echo "/$$file" | $(SED) -e 's,/[^/]*$$,,'`; \
	    if test -d "$(distdir)/$$file"; then \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -fpR $(srcdir)/$$file "$(distdir)$$dir" || exit 1; \
	      find "$(distdir)/$$file" -type d ! -perm -700 -exec chmod u+rwx {} \;; \
	    fi; \
	    cp -fpR $$d/$$file "$(distdir)$$dir" || exit 1; \
	  else \
	    test -f "$(distdir)/$$file" \
	    || cp -p $$d/$$file "$(distdir)/$$file" \
	    || exit 1; \
	  fi; \
	done

# Doc messages.
all:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) foot_message ; fi

head_message:
	@echo "=============================================================================="
	@echo " This is a directory for tests. Beside the 'make check'"
	@echo " one can use the following commands:"
	@echo "  make in:           generate all input dirs."
	@echo "  make failed-check: run check again on all directories"
	@echo "                     with missing report or failed report."
	@echo "  make complete-check: for developers, makes long and extensive tests."
	@echo "  make X.in:         generate input dir for directory X."
	@echo "  make X.check:      generate a report for directory X"
	@echo "                     (if not already existing)."
	@echo "  make X.recheck:    force the creation of the report in directory X."
	@echo "  make X.clean:      clean the given directory X."
	@echo "  make X.diff:       make the difference between output and the reference"
	@echo "                     (with the environment variable DIFF)"
	@echo "  make X.updateref   update the reference with the output"
	@echo "                     (prompt the overwrite)"

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "=============================================================================="
