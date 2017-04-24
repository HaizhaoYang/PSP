check_PROGRAMS = \
  pdgemmScaling.x \
  pzgemmScaling.x \
  testList.x \
  spBLAS.x \
  pdgemm.x \
  matSum.x \
  compareSpBLAS.x

pdgemmScaling_x_SOURCES = pdgemmScaling.F90
pdgemmScaling_x_FCFLAGS = -I../src
pdgemmScaling_x_LDADD = ../src/libpspBLAS.la

pzgemmScaling_x_SOURCES = pzgemmScaling.F90
pzgemmScaling_x_FCFLAGS = -I../src
pzgemmScaling_x_LDADD = ../src/libpspBLAS.la

compareSpBLAS_x_SOURCES = compareSpBLAS.F90
compareSpBLAS_x_FCFLAGS = -I../src
compareSpBLAS_x_LDADD = ../src/libpspBLAS.la

testList_x_SOURCES = testList.F90
testList_x_FCFLAGS = -I../src
testList_x_LDADD = ../src/libpspBLAS.la

spBLAS_x_SOURCES = spBLAS.F90
spBLAS_x_FCFLAGS = -I../src
spBLAS_x_LDADD = ../src/libpspBLAS.la

pdgemm_x_SOURCES = pdgemm.F90
pdgemm_x_FCFLAGS = -I../src
pdgemm_x_LDADD = ../src/libpspBLAS.la

matSum_x_SOURCES = matSum.F90
matSum_x_FCFLAGS = -I../src
matSum_x_LDADD = ../src/libpspBLAS.la

LOG_COMPILER = $(TESTS_CMD)

TESTS = \
  pdgemmScaling.x \
  pzgemmScaling.x \
  testList.x \
  spBLAS.x \
  pdgemm.x \
  matSum.x \
  compareSpBLAS.x

CLEANFILES = \
  zgemm.txt \
  dgemm.txt \
  zlap.txt \
  zspmm.txt \
  dspmm.txt \
  zmspm.txt \
  dlap.txt \
  dmspm.txt \
  dspmspm.txt \
  zspmspm.txt

