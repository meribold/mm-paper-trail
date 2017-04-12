C========================================================================
C	Program mcti written by Paul Beroza
C	Copyright 1991 by the Regents of the University of California.
C	All rights reserved.
C	This program is for non-commercial use.
C	Copying and distribution by author's permission only.
C	I also don't guarantee that it will work :-)
C========================================================================
C	mcti: Computes protonation vs pH (i.e., a titration curve)
C	for a protein with interacting protonatable resdiues using
C	a Monte Carlo method (Metropolis Algorithm)
C	
C	Input: 2 files) 1. intrisic pKas for the titrating sites.
C		        2. site-sites interactions in e^2/ang units.
C
C	Here is a sample shell script to run the program
C	    % mcti <<EOF
C	    % intrinsic_pKa_file        # file w/ intrinsic pKas
C	    % interactions_file         # file w/ interactions
C	    % 10000                     # num. of Monte Carlo steps
C	    % 5.0                       # starting pH
C	    % 10.0                      # final pH
C	    % 0.5                       # pH increment
C	    % min_wint  		# ! min inter for pairs (pkunits)
C	    % fract_toler  		# tolerance for reduced sites
C	    % iseed     		# if 0, then seed "randomly" 
C	    % iopt      		# 0 = full MC; 1 = red MC
C	    % log_file			# log file
C	    % EOF
C	------------------------------------------------------------
	program mcti
	implicit none
C...arrays for titrating sites:
	integer nmax           ! max number of sites
	parameter (nmax = 1000)
	integer maxsite 	! total number of sites
	real   pkint(nmax)	! intrinsic pKas of sites
	real   wint(nmax, nmax)	! interaction matrix beteen sites (e^2/ang)
C       prot(0:nmax): prot vector: either 0 or 1 for each
C		 site (0 = unprotonated, 1 = protonated)
C		 prot(0) = total protonation
	integer prot(0:nmax)
	integer resnum(nmax)	 ! residue numbers
	character*13 resnam(nmax) !residue names 
	character*1 type(nmax)   ! 'C' = cation, 'A' = anion
	real   self(nmax)		! solvation energy: f(pH,T)
	real   qunprot(nmax)	! charge of unprot state
C 	average prot. of sites
C 	aveprot(0) average tot prot.
	real   aveprot(0:nmax)
C	------------------------------------------------------
C	...reduced site monte carlo variables (r = reduced)
	real   fract_toler	! tolerance for reduced site mc.
        integer r_maxsite       ! num of reduced sites
        integer r_ptr(nmax)      ! pointer to prot() array
	integer r_resnum(nmax)
	character*13 r_resnam(nmax)
        real    r_self(nmax)         ! charge in ionzized state
        real    r_pkint(nmax)         ! charge in ionzized state
        integer r_prot(0:nmax)  ! prot. from reduced mc step
        real   r_aveprot(0:nmax)
        real    r_qunprot(nmax)         ! charge in ionzized state
        real   r_Wint(nmax,nmax)  ! W in kcal
        integer ir, jr
        logical site_fixed 
	real   w1, w2	! weights for combining errors
C	------------------------------------------------------
	character*60  interactions_outfile, log_file
C...Titration variables
	real   beta, ph		   ! beta = 1/kT	
	parameter (beta = 557.04)
	real   phlowlim, phupplim, phinc
	integer nph
	real   oldprot
	integer isite, iopt

C...Monte Carlo variables
C 	number of Monte Carlo steps. 
C  	1 step = N attempts to change a randomly selected site.
	integer mcsteps_full, mcsteps_red
C...Correlation function variables:
C...The (auto)correlation function is used to estimate the
C...sampling errors in the Monte Carlo method.
C...A sample is considered statistically independent
C...when the correlation function falls to 1/10 of it's original
C...value
	integer tau			! time var. for corr. func
	integer taumax 			! time for corr. function
	parameter (taumax = 100)
C 	c()	correlation function for each site
	real   c(0:taumax,0:nmax) 	! cor function
	real   s(0:taumax,0:nmax)  ! array for cor. function calc.
	real   save(0:nmax) 	! array for cor. function calc.
	real   error 		! std. dev. for ave prot.
	real   prot_error(0:nmax)
	integer iold(0:taumax,0:nmax) ! array of old values

C...Energy variables:
	real   energy, eave 		! energy and average energy
	real   cenergy(0:taumax)		! correlation function for energy
	real   seave, se(0:taumax) 	! array for corr. function
	real   eold(0:taumax)		! array for corr. function
	real   grounde		! energy with all sites unprotonated

C...Two Site Transition Variables
C...Pairs of strongly interacting sites are allowed to simultaneously
C...change there protonation states.
	integer maxpairs		! max. num of pairs
	parameter (maxpairs = 500)
	integer npairs  	   ! number of strongly coupled sites
	integer pairs(2,maxpairs)  ! site numbers of coupled sites
	real   min_wint		! min inter (pK units) for 2-site transitions
        integer r_maxpairs, r_npairs
        parameter (r_maxpairs = maxpairs)
        integer r_pairs(2,r_maxpairs)
C...Misc variables:
	integer i, j
	character*60	cc1, cc2
	integer iseed	! rn gen. seed
	integer time	! returns time in integers since the epoch
			! used to seed random number generator
	integer t0, t1, t2, mclock	! timing variables
	integer istep, nstep ! GMU variables for pH loop
C========================================================================
C========================================================================
	cc1='mcti written by Paul Beroza'
	cc2='Copyright 1991 by the Regents of the Univ. of California.'
	t0 = mclock()
C	------------------------------------------------------------
C	...get the interactions and intrinsic pKas
	call getdata(nmax, wint, pkint, qunprot, maxsite,
     &	             resnum, resnam, type)
C	------------------------------------------------------------
C	...read the titration variables
	read(unit=5, fmt=*) mcsteps_full
	read(unit=5, fmt=*) mcsteps_red
	read(unit=5, fmt=*) phlowlim
	read(unit=5, fmt=*) phupplim
	read(unit=5, fmt=*) phinc
	read(5,*) min_wint ! min inter for pairs (pkunits)
	read(5,*) fract_toler
	read(5,*) iseed   ! if 0, then seed "randomly"
	read(5,*) iopt    ! 0 = full MC; 1 = red MC; 2 = Cluster
	read(5,'(a)') log_file
	open(unit=12,file=log_file,status='unknown')
C	------------------------------------------------------------
	if (iopt .gt. 2 .or. iopt .lt. 0) stop 'bad iopt'
c	...seed random number generator
	if (iseed .eq. 0) then
C Never use 0 in f77 it does not work properly
C	   iseed = time(%VAL(0))    !%VAL because time is c
           write(12,*) 'SEED not set -Dont do this!!!!'
           iseed=45543534 
	endif
        write(12,*) 'SEED: ', iseed
        call srand(iseed)              ! seed random number
C	------------------------------------------------------------
C	Full Monte Carlo
C	...calculate the ground energy of the system (no protons)
	grounde = 0
	do j = 1, maxsite
	 do i = 1, maxsite
	     grounde = grounde + (qunprot(i) * qunprot(j)) * wint(i,j)
	 enddo
	enddo
	grounde = 0.5*grounde
C	write(6,*)'grounde = ', grounde
C	------------------------------------------------------------
C	...Generate pair list for 2-site transitions.
	call get_pairs(nmax, maxsite, Wint, min_wint, beta,
     &	               maxpairs, npairs, pairs)
	write(12,*) 'Full MC. pairs: Wint <',  min_wint
	do i = 1, npairs
	   write(12,500) resnam(pairs(1,i)), resnum(pairs(1,i)), 
     &	              resnam(pairs(2,i)), resnum(pairs(2,i)), 
     &		      beta*wint(pairs(1,i), pairs(2,i))/log(10.0)
	enddo
500	format(a13,1x,i4,' <-> ',a13,1x,i4,' = ',f6.1)
C	------------------------------------------------------------
C	------------------------------------------------------------
C	...pH LOOP
	nph = (phupplim-phlowlim)/phinc +1
c gmu	nph = 0
c gmu	do ph = phlowlim, phupplim+0.001, phinc  ! PH:START
c gmu	   nph = nph + 1
c gmu	enddo
	write(6,*) maxsite, nph
	do i = 1, maxsite
	   write(6,'(f10.5,1x,a1,1x,a13,2x,i4)')
     &			pkint(i), type(i), resnam(i), resnum(i)
	enddo
c gmu	DO ph = phlowlim, phupplim, phinc  ! PH:START
	DO istep = 0, nph, 1
	   ph = phlowlim + istep*phinc
	   t1 = mclock()
C	  ----------------------------------------------------------
	   call init(nmax, prot, maxsite)   ! ...initialize prot state
C	   call showlat(nmax, prot, maxsite)
C	  ----------------------------------------------------------
	   call mc(nmax, prot, wint, pkint, self, qunprot, beta, ph,
     &	         maxsite, mcsteps_full, aveprot, taumax, c, s, 
     &           save, iold, energy, eave, cenergy, seave, se,
     &		 eold, maxpairs, npairs, pairs)
	   t2 = mclock()
	   write(12,'("TIME: Full Monte Carlo (min): ",2(i8,1x),f10.2)')
     &	          maxsite, mcsteps_full, float(t2 - t1)/(60.0*100.0)
C	  ----------------------------------------------------------
	   write(unit=6, fmt=*) ph
	   write(unit=12, fmt=*) 'pH: ', ph
C	   -------------------------------------------------------
c          ...calc. the std dev. for each prot. point
	   if (iopt .ne. 0) write(12,*) 'Prot after full Monte Carlo'
	   do i = 0, maxsite
	    call geterr(taumax, c(0,i), mcsteps_full,error)
	    prot_error(i) = error
	    if (iopt .eq. 0) then
	       write(unit=6, fmt=*) i, aveprot(i), prot_error(i)
	    else
	       write(unit=12, fmt=*) i, aveprot(i), prot_error(i)
	    endif
	   enddo
C	   -------------------------------------------------------
C	   ...Reduced Site Techniques
	   if (iopt .ne. 0) then 
C	      ...get reduced set
	      call reduce(fract_toler, beta, nmax,
     &                    maxsite, wint, pkint, resnum,
     &                    resnam, qunprot, aveprot,
     &                    r_maxsite, r_Wint, r_pkint, r_resnum,
     &                    r_resnam, r_qunprot, r_ptr)
C	      ------------------------------------------------------------
	      write(12,*) 'Reduced sites at pH:', ph
	      do i = 1, r_maxsite
		 write(12,'(a13,1x,i4,1x,f6.1)') 
     &		         r_resnam(i), r_resnum(i), r_pkint(i)
	      enddo
	      if (iopt .eq. 1) then
C	      =========================================================
C	      ...Reduced Monte Carlo
C	         ------------------------------------------------------------
C	         ...Generate pair list for 2-site transitions.
	         call get_pairs(nmax, r_maxsite, r_Wint, min_wint, beta,
     &	               r_maxpairs, r_npairs, r_pairs)
	         write(12,*) 'Reduced pairs pH:', pH
	         do i = 1, r_npairs
	            write(12,500) r_resnam(r_pairs(1,i)), 
     &		      r_resnum(r_pairs(1,i)), r_resnam(r_pairs(2,i)), 
     &			r_resnum(r_pairs(2,i)), 
     &		      beta*r_Wint(r_pairs(1,i), r_pairs(2,i))/log(10.0)
	         enddo
C	         ------------------------------------------------------------
	         call init(nmax, r_prot, r_maxsite)   ! ...initialize prot state
C	         ------------------------------------------------------------
C	         ...note: for now the energy stuff has no meaning
		 t1 = mclock()
	         call mc(nmax, r_prot, r_Wint, r_pkint, r_self, r_qunprot,
     &	              beta, ph, r_maxsite, mcsteps_red, r_aveprot, 
     &		      taumax, c, s, 
     &                save, iold, energy, eave, cenergy, seave, se,
     &			eold, r_maxpairs, r_npairs, r_pairs)
		 t2 = mclock()
		 write(12,'("TIME: Red. M.C. (min): ",2(i8,1x),f10.2)')
     &               r_maxsite, mcsteps_red, float(t2 - t1)/(60.0*100.0)

	         do i = 1, r_maxsite  ! exclude total prot.
	            isite = r_ptr(i)
	            call geterr(taumax, c(0,i), mcsteps_red, error)
		    if (error .lt. 1.0e-8) then
		       aveprot(isite) = r_aveprot(i)
		       prot_error(isite) = error
		    else
		       w1 = 1.0d0/(prot_error(isite)*prot_error(isite))
	               w2 = 1.0d0/(error*error)
		       aveprot(isite) =
     &		         (w1*aveprot(isite) + w2*r_aveprot(i))/(w1 + w2)
	               prot_error(isite) = sqrt(1.0d0/(w1 + w2))
		    endif
	         enddo
	         aveprot(0) = 0.0d0
	         prot_error(0) = 0.0d0
	         do i = 1, maxsite
	            aveprot(0) =   aveprot(0) + aveprot(i)
	            prot_error(0) = prot_error(0) + 
     &	                      prot_error(i)*prot_error(i)
	         enddo
	         prot_error(0) = sqrt(prot_error(0))
	         do i = 0, maxsite
	            write(unit=6, fmt=501) i, aveprot(i), prot_error(i)
 501		    format (i4, f10.5, f10.5)
	         enddo
	      else 
	         stop ' bad IOPT'
	      endif
	   endif  !IOPT .ne. 0
C=========================================================
	ENDDO ! PH END
	t1 = mclock()
	write(12,'("TIME: Total (min): ",f10.2)')
     &                   float(t1 - t0)/(60.0*100.0)
	end
C========================================================================
C------------------------------------------------------------------------
C subroutine geterr(taumax,cc,mcsteps,error): computes the standard
C	deviation from the correlation function.  If the correlation
C	does not fall to 1/10 the value of c(0) then assign 999 to the
C	error.
C 
	subroutine geterr(taumax,cc,mcsteps,error)
	implicit none
	integer taumax
	real   cc(0:taumax)   	! correlation function for 1 site
	integer mcsteps     	! num of Monte Carlo steps
	real   error		! error (returned)
	logical flag
	real   ctol
	integer i, ict

C...flag = .true. when the  correlation function is small
	flag = .false.
	error = 0
	ctol = 0.1*cc(0)

C 
C... if cc(0) is aleady small, the correlation is = 0
	if (abs(cc(0)) .lt. 1.0e-5) then
	    ict = 0
	    flag  = .true.
	endif

	do 100 i = 1, taumax
		if (cc(i) .lt. ctol .and. .not. flag ) then
		   ict = i
		   flag = .true.
		endif
100	continue
	if ( .not. flag ) then   ! correlations too strong
		error = 999
	else
C.... mcsteps/ict  = num. of statistically independent samples
		error = sqrt( (cc(0)*ict) / float(mcsteps) )
	endif
	return
	end

