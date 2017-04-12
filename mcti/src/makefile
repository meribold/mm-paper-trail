FC = gfortran -O3
OBJS = mcti.o mc.o mcsubs.o step.o pair_flip.o get_pairs.o reduce.o
SRCS = mcti.f mc.f mcsubs.f step.f pair_flip.f get_pairs.f reduce.f

mcti: $(OBJS)
	$(FC) -o mcti $(OBJS)

