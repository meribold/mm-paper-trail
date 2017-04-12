C========================================================================
C	Program mcti written by Paul Beroza
C	Copyright 1991 by the Regents of the University of California.
C	All rights reserved.
C	This program is for non-commercial use.
C	Copying and distribution by author's permission only.
C	I also don't guarantee that it will work :-)
C========================================================================
C subroutine getdata:  reads in the intrinsic pKas and site-site
C			interactions.
      subroutine getdata(nmax, wint, self, qunprot, maxsite,
     *      resnum, resnam, type)
      implicit none
      integer nmax, maxsite
      real   wint(nmax, nmax), self(nmax), qunprot(nmax)
      integer resnum(nmax)
      character*13 resnam(nmax)
      character*1 type(nmax)
      integer i, j, icheck, jcheck
      character pkfile*30, gfile*30
      maxsite = 0
      do 5 i = 1, nmax
         self(i) = 0
         qunprot(i) = 0
         do 5 j = 1, nmax
            wint(i,j) = 0.
    5 continue
      read(5 ,'(a)') pkfile
      read(5, '(a)') gfile
      open(unit=1, file=pkfile, status='old') 
      open(unit=2, file=gfile, status='old') 
      do 10 i = 1, nmax
	 read(unit=1, fmt=*, end=20) self(i), type(i), resnam(i) 
         resnum(i)=i
         if (type(i) .eq. 'A') then
            qunprot(i) = -1
         else
            qunprot(i) = 0
         end if
	 goto 10
   10 continue
      stop 'nmax exceeded'
   20 continue
      maxsite = i - 1
C      write(6,*) maxsite
C      do i = 1, maxsite
C         write(6,500) self(i), type(i), resnam(i), resnum(i)
C      enddo
500   format (f10.5,1x,a1,1x,a3,2x,i4)
      do 30 i = 1, maxsite
         do 30 j = 1, maxsite
            read(unit=2, fmt=*) icheck, jcheck, wint(i,j)
            if ((i .ne. icheck) .or. (j .ne. jcheck)) then
               stop 'number mismatch in g.out file'
            end if
   30 continue
      close(unit=1) 
      close(unit=2) 
      return 
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c.. subrountine init : initializes protonation to a random state
      subroutine init(nmax, lattice, maxsite)
      implicit none
      integer nmax, maxsite
      integer lattice(0:nmax)
      integer i
      real    rand

      lattice(0) = 0
      do 10 i = 1, maxsite
         if (rand() .gt. 0.5d0) then
            lattice(i) = 1
         else
            lattice(i) = 0
         end if
         lattice(0) = lattice(0) + lattice(i)
   10 continue
      return 
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C subroutine showlat:  writes out lattice (for debugging)
      subroutine showlat(nmax, lattice, maxsite)
      implicit none
      integer nmax, maxsite, lattice(0:nmax)
      integer i

      write(12,'(i4,1x,200i1)') (lattice(i),i=0,maxsite)
      return 
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C subroutine showq: writes out charges (for debugging)
      subroutine showq(nmax, prot, qunprot, maxsite)
      implicit none
      integer nmax, maxsite, prot(0:nmax)
      real   qunprot(nmax)
      integer i
      do 10 i = 1, maxsite
         write(unit=6, fmt='(i3,1x,f4.1)') i, qunprot(i) + prot(i)
   10 continue
      return 
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C subroutine fillself: calculates energy to protonate a site, in the
C			absence of interactions with other sites.
      subroutine fillself(beta, ph, pkint, nmax, maxsite, self)
      implicit none
      integer nmax, maxsite
      real   beta, ph, pkint(nmax), self(nmax)
      integer i
      real   const
      const = - (log(10.0) / beta)
      do 10 i = 1, maxsite
         self(i) = const * (pkint(i) - ph)
   10 continue
      return 
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C function deltae: calculates change in energy upon changes the
C		   protonation at site iflip
      function deltae(nmax, maxsite, self, prot, qunprot, wint, iflip)
      implicit none
      integer nmax, maxsite, prot(0:nmax), iflip
      real   self(nmax), qunprot(nmax), wint(nmax, nmax), deltae
      integer dp, j

      dp = 1 - (2 * prot(iflip))
      deltae = 0
      do 10 j = 1, maxsite
         deltae = deltae + ((dp * (qunprot(j) + prot(j))) * 
     &			wint(iflip,j))
   10 continue
      deltae = deltae + (dp * self(iflip))
      return 
      end
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C subroutine thermal: mcsteps of Monte Carlo to approach
C			equilibium before data is taken.
C	Note: This currently does not allow a 2 site transition.
      subroutine thermal(nmax, maxsite, self, prot, qunprot, wint, beta
     &, mcsteps)
      implicit none
      integer nmax
c prot vector: either up or down
      integer prot(0:nmax)
c interaction matrix beteen sites
      real   wint(nmax, nmax)
c energy of solvation = f(pH,T)
      real   self(nmax)
c array of charges of unprot state
      real   qunprot(nmax)
c number of cycles to run 
      integer mcsteps
c beta = 1/kT
      real   beta
c total number of sites
      integer maxsite
      integer i, iflip
c	   ...loop mcsteps (mcsteps*maxsite) to thermalize
c delta e from given flip
      real   de, deltae
      real   rand
      do 90 i = 1, mcsteps * maxsite
c choose spin
         iflip = int((maxsite * rand()) + 1)
c		...decide
         de = deltae(nmax,maxsite,self,prot,qunprot,wint,iflip)
         if (de .lt. 0) then
            prot(iflip) = 1 - prot(iflip)
            prot(0) = (prot(0) + (2 * prot(iflip))) - 1
         else if (exp(- (beta * de)) .gt. rand()) then
            prot(iflip) = 1 - prot(iflip)
            prot(0) = (prot(0) + (2 * prot(iflip))) - 1
         else
         end if
   90 continue
      return 
      end

