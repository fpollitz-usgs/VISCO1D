# Makefile for DECAY, DECAY4, DECAY4M, VSPHDEP,
# VSPHM, VTORDEP, and STRAINA
# Options below are for Absoft compiler on MAC powerbook OS10.2, OS10.3
#FC=f77
#FFLAGS=-O -q -W -s
# Options below are for ifort compiler on Intel IMAC OSX10.4.11
#FC=ifort
#FFLAGS=-O2 -ftz -ip -ipo -axP -mp -align all -extend_source
# Options below are for GNU compiler
FC=gfortran
FFLAGS=-ffixed-form -ffixed-line-length-132 

all : decay decay4 decay4m vsphdep vsphm vtordep strainA

decay : decay.f
	${FC} -o $@ ${FFLAGS} decay.f

decay4 : decay4.f
	${FC} -o $@ ${FFLAGS} decay4.f

decay4m : decay4m.f bsubs.f
	${FC} -o $@ ${FFLAGS} decay4m.f bsubs.f

strainA : strainA.f lgndrl.f splneq.f
	${FC} -o $@ ${FFLAGS} strainA.f lgndrl.f splneq.f

vsphdep : vsphdep.f
	${FC} -o $@ ${FFLAGS} vsphdep.f

vsphm : vsphm.f bsubs.f bsubs2.f splneq.f
	${FC} -o $@ ${FFLAGS} vsphm.f bsubs.f bsubs2.f splneq.f

vtordep : vtordep.f
	${FC} -o $@ ${FFLAGS} vtordep.f

