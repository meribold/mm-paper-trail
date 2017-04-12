C========================================================================
C	Program mcti written by Paul Beroza
C	Copyright 1991 by the Regents of the University of California.
C	All rights reserved.
C	This program is for non-commercial use.
C	Copying and distribution by author's permission only.
C	I also don't guarantee that it will work :-)
C========================================================================
c	subroutine step  take 1 monte carlo step
c		... involves maxsite attempts to change prot
      subroutine step(nmax, maxsite, self, prot, qunprot, wint, 
     &      		beta, energy, maxpairs, npairs, pairs)
      implicit none
c max number of sites
      integer nmax
c total number of sites
      integer maxsite
c energy of solvation = f(pH,T)
      real self(nmax)
c prot(0) = total protonation
c prot vector: either up or down
      integer prot(0:nmax)
c array of charges of unprot state
      real qunprot(nmax)
c interaction matrix beteen sites
      real wint(nmax, nmax)
c beta = 1/kT
      real beta
      real energy
      integer i, iflip
      real de, deltae
      real rand
C
C...2site transition variables
	integer maxpairs
	integer npairs
	integer pairs(2, maxpairs)

C... in this sampling we are considering the pairs of strongly interacting
C... residues as separate "sites" that can undergo change (i.e., a 2site
C... transition), 
      do 100 i = 1, maxsite + npairs
c choose spin
         iflip = int(((maxsite + npairs)* rand()) + 1)
c		...decide
	 if (iflip .gt. maxsite) then  ! we will attempt a 2-site transition
	    iflip = iflip - maxsite
	    call pair_flip(pairs(1,iflip), energy,
     &	    beta, nmax, maxsite, prot, qunprot, wint, self)
	 else	! a 1-site transition
            de = deltae(nmax,maxsite,self,prot,qunprot,wint,iflip)
C	    ...now the standard Metropolis criteria:
C		if de < 0 change
C		if de >= 0 change with probability exp(-beta*de)
            if (de .lt. 0) then
               prot(iflip) = 1 - prot(iflip)
c gmu               prot(0) = (prot(0) + (2 * prot(iflip))) - 1
               energy = energy + de
            else if (exp(- (beta * de)) .gt. rand()) then
               prot(iflip) = 1 - prot(iflip)
c gmu              prot(0) = (prot(0) + (2 * prot(iflip))) - 1
               energy = energy + de
            end if
	 endif
  100 continue
      prot(0) = 0 
      do i = 1, maxsite 
         prot(0) = prot(0) + prot(i)
      enddo
      return 
      end

