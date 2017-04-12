C========================================================================
C       Program mcti written by Paul Beroza
C       Copyright 1991 by the Regents of the University of California.
C       All rights reserved.
C       This program is for non-commercial use.
C       Copying and distribution by author's permission only.
C       I also don't guarantee that it will work :-)
C========================================================================
C       -----------------------------------------------------------------
        subroutine get_pairs(nmax, maxsite, Wint, min_wint, beta,
     &                 maxpairs, npairs, pairs)
        implicit none
        integer nmax, maxsite
        real Wint(nmax, nmax), min_wint, beta
        integer maxpairs, npairs, pairs(2,maxpairs)
        integer i, j
	real min_g

	min_g = min_wint * log(10.0) / beta
C       write(6,*) 'min_wint = ', min_wint
        npairs = 0
        do i = 1, maxsite
           do j = i + 1, maxsite
              if (wint(i,j) .gt. min_g) then
                 npairs =   npairs + 1
                 if (npairs .gt. maxpairs) stop ' too many pairs'
                 if (npairs + maxsite .gt. nmax) stop ' nmax too small'
                 pairs(1,npairs) = i
                 pairs(2,npairs) = j
              endif
           enddo
        enddo
	return
	end
C       -----------------------------------------------------------------

