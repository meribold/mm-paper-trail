C========================================================================
C	Program mcti written by Paul Beroza
C	Copyright 1991 by the Regents of the University of California.
C	All rights reserved.
C	This program is for non-commercial use.
C	Copying and distribution by author's permission only.
C	I also don't guarantee that it will work :-)
C========================================================================
C subroutine reduce:  
C	 	  reduces the number of sites in the Monte Carlo
C		  calculation to only those which have some
C		  fractional protonation
C		  
C		  
	subroutine reduce(fract_toler, beta, nmax, 
     &	   		  maxsite, wint, pkint, resnum,
     &			  resnam, qunprot, aveprot,
     &	   		  r_maxsite, r_Wint, r_pkint, r_resnum,
     &			  r_resnam, r_qunprot, r_ptr)
	implicit none
	real    fract_toler		! fract. prot. tolerance
	real   beta
	integer nmax 			! max number of sites
C	... full site variables
	integer maxsite 	! total number of sites
	real   wint(nmax, nmax) ! interaction matrix beteen sites
	real   pkint(nmax)      ! pKintr energy, 
	integer resnum(nmax)
	character*13 resnam(nmax)
	real   qunprot(nmax) 	! array of charges of unprot state
	real   aveprot(0:nmax) 	! aveprot(0) tot prot. ave
C	--------------------------------------------
C	...arrays for second set of statistics (r = reduced)
	integer r_maxsite
	integer r_ptr(nmax)      ! pointer to prot() array
	real    r_qunprot(nmax)         ! charge in ionzized state
	real   r_Wint(nmax,nmax)  ! W in kcal
	real   r_pkint(nmax)
	integer r_resnum(nmax)
	character*13 r_resnam(nmax)
C	-------------------------------------------
C	...local vars
	integer ir, jr, j, i, isite, jsite
	logical site_fixed 
C	===================================================
C	--------------------------------------------------
C	... call mc again with reduced set
C	write(12,*) 'setting up arrays for reduced mc' 
	r_maxsite = 0
C	----------------------------------------
	do j = 1, maxsite
	   if (aveprot(j) .gt. 1.0d0 - fract_toler .or. 
     &	       aveprot(j) .lt. fract_toler) then
	   else
	      r_maxsite = r_maxsite + 1
	      if (r_maxsite .gt. nmax) stop ' r_maxsite .gt. nmax'
	      r_ptr(r_maxsite) = j
	      r_qunprot(r_maxsite) = qunprot(j)
	      r_pkint(r_maxsite) = pkint(j)
	      r_resnum(r_maxsite) = resnum(j)
	      r_resnam(r_maxsite) = resnam(j)
	   endif
	enddo
	WRITE(12,*) 'N REDUCED SITES: ', r_maxsite
C	------------------------------------------------ 
C	...fill W
	do j = 1, r_maxsite
	   jsite = r_ptr(j)
	   do i = j, r_maxsite
	      isite = r_ptr(i)
	      r_Wint(i,j) = Wint(isite,jsite) 
	      r_Wint(j,i) = r_Wint(i,j)
	   enddo
	enddo
C	------------------------------------------------ 
C	...adjust r_pkint for fixed charges
	do ir = 1, r_maxsite
	   isite = r_ptr(ir)
	   do jsite = 1, maxsite
	      site_fixed = .true.
	      do jr = 1, r_maxsite
		 if (r_ptr(jr) .eq. jsite) site_fixed = .false.
	      enddo
	      if (site_fixed) 
     &		      r_pkint(ir) = r_pkint(ir) 
     &			   - (qunprot(jsite) + aveprot(jsite))*
     &			     ( (beta*Wint(isite,jsite))/log(10.0d0) )
	   enddo
	enddo
C	------------------------------------------------ 
	return
	end

