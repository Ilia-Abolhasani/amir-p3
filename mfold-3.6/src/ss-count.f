c     This is a simplified redo of John Jaeger's sscount.f
c     A series of RNA/DNA structures in .ct format are read
c     via standard input.
c     It is assumed that all the structures are on the same sequence.
c     Standard output: The first record contains the number of foldings.
c     Record i+1 contains: hstnum(i) j, where j is the number of times 
c     the ith base is single stranded in all the  foldings.

      parameter (maxsiz=20000)
      integer count(maxsiz),nbases,hstnum(maxsiz),check
      character*1 seq(maxsiz)

c     k counts the nummber of foldings
      k = 0
      do while (1.eq.1)
         k = k + 1
         if (k.eq.1) then
            read(5,*,end=99) nbases
            do n = 1,nbases
               count(n) = 0
            enddo
         else
            read(5,*,end=98) check
            if (check.ne.nbases) then
               write(6,*) 'STOP: Incorrect header in ct file.'
               call exit(1)
            endif
         endif
         do i = 1,nbases
            read(5,*,end=97) check,seq(i),itmp,itmp,j,hstnum(i)
            if (check.ne.i) then
               write(6,*) 'STOP: corrupted ct file.'
               call exit(1)
            endif
            if (j.eq.0) count(i) = count(i) + 1
         enddo
      enddo

 99   write(6,*) 'STOP: No structures.'
      call exit(1)
 97   write(6,*) 'STOP: Truncated .ct file.'
      call exit(1)
 98   k = k - 1
      write(6,*) k
      do i = 1,nbases
         write(6,*) hstnum(i),count(i),' ',seq(i)
      enddo
      call exit(0)
      end
