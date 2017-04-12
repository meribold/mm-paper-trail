C========================================================================
C	Program mcti written by Paul Beroza
C	Copyright 1991 by the Regents of the University of California.
C	All rights reserved.
C	This program is for non-commercial use.
C	Copying and distribution by author's permission only.
C	I also don't guarantee that it will work :-)
C========================================================================
C subroutine mc:  carries out the monte carlo sampling for a given
C	 	  pH value and returns the average values of the 
C		  protonation for each site and the correlation functions
C		  used to compute the energy.
	subroutine mc(nmax, prot, wint, pkint, self, qunprot, beta, ph, 
     &    maxsite, mcsteps, aveprot, taumax, c, s, save, iold, 
     &	  energy, eave, ce, seave, se, eold, maxpairs, npairs, pairs)
	implicit none
	integer nmax 			! max number of sites
	integer prot(0:nmax)    ! prot vector: either up or down
			        ! prot(0) = total protonation
	real   wint(nmax, nmax) ! interaction matrix beteen sites
	real   pkint(nmax)      ! pKintr energy, 
	real   self(nmax) 	! energy of solvation = f(pH,T)
	real   qunprot(nmax) 	! array of charges of unprot state
	real   beta, ph 	! beta = 1/kT
	integer maxsite 	! total number of sites
	integer mcsteps 	! num of mc steps
	real   aveprot(0:nmax) 	! aveprot(0) tot prot. ave
				! array of prot. ave for sites
	integer taumax
	real   c(0:taumax,0:nmax) 	! cor function
	real   save(0:nmax) 		! array for cor. function calc.
	real   s(0:taumax,0:nmax) 	! array for cor. function calc.
	integer iold(0:taumax,0:nmax)  ! array of old values
	integer t, i, j 			! t is mc carlo unit of time
C	--------------------------------------------
c	...vars for energy corr function and ave.
	integer tau 	! time var. for corr. func
	real   energy, eave
	real   ce(0:taumax)
	real   seave, se(0:taumax)
	real   eold(0:taumax)
	real   esum
C	--------------------------------------------
C	...2site transition variables:
	integer maxpairs
	integer npairs
	integer pairs(2, maxpairs)
C	===================================================

C	--------------------------------------------------
	eave = 0.0
	seave = 0.0
	do tau = 0, taumax
	   ce(tau) = 0.0
	   se(tau) = 0.0
	   eold(tau) = 0.0
	enddo
	do j = 0, maxsite
	   prot(j) = 0.0
	   aveprot(j) = 0.0
	   save(j) = 0.0
	   do tau = 0, taumax
	      c(tau,j) = 0.0
	      s(tau,j) = 0.0
	      iold(tau,j) = 0.0
	   enddo
	enddo
C	--------------------------------------------------
C	...compute the self energy to protonate each site
	call fillself(beta, ph, pkint, nmax, maxsite, self)
c 	...thermalization:  do mcsteps of monte carlo before sampling
	call thermal(nmax, maxsite, self, prot, qunprot, wint, beta, 
     &			mcsteps)
	energy = esum(nmax,maxsite,self,prot,qunprot,wint)
C	--------------------------------------------------------
c	...statistics loop: take 1 mc step and calcalate averages
C	call showlat(nmax, prot, maxsite)
	do t = 1, mcsteps
C	   ------------------------------------------------------
	   call step(nmax, maxsite, self, prot, qunprot, 
     &		   wint, beta, energy, maxpairs, npairs, pairs)
C	   call showlat(nmax, prot, maxsite)
	   eave = eave + energy
	   seave = seave + energy
	   se(0) = se(0) + (energy * energy)
	   do tau = 1, taumax
	      se(tau) = se(tau) + (energy * eold(tau))
	   enddo
	   do tau = taumax, 2, -1
	      eold(tau) = eold(tau - 1)
	   enddo
	   eold(1) = energy
C	   ------------------------------------------------------
	   do j = 0, maxsite
	      aveprot(j) = aveprot(j) + prot(j)
	      IF (aveprot(j) .lt. 0.0d0) 
     &	         write(6,*)'aveprot neg. aveprot, prot', 
     &	            aveprot(j), prot(j),'---',j,t,maxsite, ph
	      save(j) = save(j) + prot(j)
	      s(0,j) = s(0,j) + (prot(j) * prot(j))
	      do tau = 1, taumax
	        s(tau,j) = s(tau,j) + (prot(j) * iold(tau,j))
	      enddo
	      se(tau) = se(tau) + (energy * eold(tau))
	      do tau = taumax, 2, -1
	         iold(tau,j) = iold(tau-1,j)
	      enddo
	      iold(1,j) = prot(j)
	   enddo
C	   ------------------------------------------------------
	enddo
C	--------------------------------------------------------
	eave = eave / float(mcsteps)
	seave = (seave/mcsteps) * (seave/mcsteps)
	do tau = 0, taumax
	   ce(tau) = ((1.0 / float(mcsteps - tau)) * se(tau)) - seave
	enddo
C	--------------------------------------------------------
	do j = 0, maxsite
	   aveprot(j) = aveprot(j) / float(mcsteps)
	   save(j) = (save(j)/mcsteps) * (save(j)/mcsteps)
	   do tau = 0, taumax
	      c(tau,j) = ((1.0 / float(mcsteps - tau)) * s(tau,j))
     &			 - save(j)
	   enddo
	enddo
C	--------------------------------------------------------
	return 
	end
C----------------------------------------------------------------------
	function esum(nmax, maxsite, self, prot, qunprot, wint)
	implicit none
	real   esum
	integer nmax, maxsite, prot(0:nmax)
	real   self(nmax), qunprot(nmax), wint(nmax, nmax)
	integer i, j
	external showlat	
 
C	call showlat(nmax, prot, maxsite)
	esum = 0
	do j = 1, maxsite
	   esum = esum + (prot(j) * self(j))
	   do i = 1, maxsite
	      esum = esum + 
     &	        0.5*(
     &	              ((qunprot(i) + prot(i)) * (qunprot(j) + prot(j)))
     &	                                -
     &                (qunprot(i) * qunprot(j))
     &	             ) * wint(i,j)
	   enddo
	enddo
C	write(6,*) 'initial energy = ', esum
	return 
	end

