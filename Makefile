# ============================================================================
#  SU3HOB-ladder -- harmonic-oscillator (Talmi-Moshinsky) brackets in the
#  SU(3) basis, with the SU(3) > SO(3) isofactors generated internally from
#  the U(2) pseudo-spin ladder operators.  Self-contained: no external SU(3)
#  coupling library is required.
#
#    make            build the validation driver and the CLI
#    make check      run the validation driver against test_hob_su3.ref
#    make bench      build the timing / accuracy drivers in bench/
#    make clean
#
#  BLAS/LAPACK is needed only by the bench/ drivers (dsyev, matmul).
# ============================================================================
# Make predefines FC = f77, so `?=' would never take effect; replace the
# built-in default only, leaving a user-supplied FC (command line or
# environment) alone.
ifeq ($(origin FC),default)
FC       = gfortran
endif
FFLAGS  ?= -O2 -Wall -Wno-unused-dummy-argument
LAPACK  ?= -llapack -lblas

OBJ  = WignerSymbol.o hob_su3.o
BENCH = bench/ladder

all: test_hob_su3 hob

test_hob_su3: $(OBJ) tmb_kam.o test_hob_su3.o
	$(FC) $(FFLAGS) -o $@ $(OBJ) tmb_kam.o test_hob_su3.o

hob: $(OBJ) hob.o
	$(FC) $(FFLAGS) -o $@ $(OBJ) hob.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

hob_su3.o: WignerSymbol.o
test_hob_su3.o: hob_su3.o tmb_kam.o
hob.o: hob_su3.o

# The driver also prints a wall-clock section, which cannot match across
# machines; the reference holds the deterministic part only, so compare that.
check: test_hob_su3
	@./test_hob_su3 > test_hob_su3.out
	@sed '/^Timing/,$$d' test_hob_su3.out > test_hob_su3.chk
	@if diff -q test_hob_su3.chk test_hob_su3.ref >/dev/null 2>&1; then \
	   echo "PASS: numerical output matches test_hob_su3.ref"; \
	 else \
	   echo "FAIL: differs from reference"; diff test_hob_su3.ref test_hob_su3.chk | head; \
	 fi

bench: $(BENCH)

bench/ladder: $(OBJ) bench/ladder.f90
	$(FC) $(FFLAGS) -o $@ bench/ladder.f90 $(OBJ) $(LAPACK)

clean:
	rm -f *.o *.mod test_hob_su3 hob $(BENCH) test_hob_su3.out test_hob_su3.chk

.PHONY: all check bench clean
