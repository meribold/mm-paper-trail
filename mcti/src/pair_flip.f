C========================================================================
C	Program mcti written by Paul Beroza
C	Copyright 1991 by the Regents of the University of California.
C	All rights reserved.
C	This program is for non-commercial use.
C	Copying and distribution by author's permission only.
C	I also don't guarantee that it will work :-)
C========================================================================
C subroutine pair_flip: attempt a simultaneous change in the protonation
C			state of 2 sites, and accept the change based on
C			the standard Metropolis criteria.
	subroutine pair_flip(pair, energy,
     &	   beta, nmax, maxsite, prot, qunprot, wint, self)
	implicit none
	real   energy
	integer pair(2)		!sites to change
	real   beta		! 1/kT
	integer nmax, maxsite
	integer prot(0:nmax)	! protonation state
	real   qunprot(nmax)	! charged state
	real   wint(nmax, nmax)   ! e. s. interaction
	real   self(nmax)		! intrinsic E to protonate site

	integer newprot(2)	! new prot state of sites
	integer dp(2)		! change in prot state of sites
	integer i
	real   de, rand

	do i = 1, 2
           if (rand() .gt. 0.5d0) then
              newprot(i) = 1
           else
              newprot(i) = 0
           end if
c gmu	   newprot(i) = int(2*rand())
c gmu	   if (newprot(i) .gt. 1) stop ' bad pair_flip'
	   dp(i) = newprot(i) - prot(pair(i))
	enddo

	if (dp(1) .eq. 0 .and. dp(2) .eq. 0) return

	de = dp(1)*self(pair(1)) + dp(2)*self(pair(2))

C...change in self energy
	do i = 1, maxsite
	   if (i .ne. pair(1) .and. i .ne. pair(2)) then
	      de = de + 
     &	      dp(1) * (qunprot(i) + prot(i)) * wint(i,pair(1))
     &		      +
     &	      dp(2) * (qunprot(i) + prot(i)) * wint(i,pair(2))
	   endif
	enddo

C...change in interaction energy
	de = de + wint(pair(1), pair(2)) * (
     &           (qunprot(pair(1)) + newprot(1))*
     &           (qunprot(pair(2)) + newprot(2)) -
     &           (qunprot(pair(1)) + prot(pair(1)))*
     &           (qunprot(pair(2)) + prot(pair(2))) )

	if (de .lt. 0) then
	   prot(pair(1)) = prot(pair(1)) + dp(1)
	   prot(pair(2)) = prot(pair(2)) + dp(2)
	   prot(0) = prot(0) + dp(1) + dp(2)
	   energy = energy + de
	elseif (exp(-beta*de) .gt. rand()) then
	   prot(pair(1)) = prot(pair(1)) + dp(1)
	   prot(pair(2)) = prot(pair(2)) + dp(2)
	   prot(0) = prot(0) + dp(1) + dp(2)
	   energy = energy + de
	endif
	return
	end

