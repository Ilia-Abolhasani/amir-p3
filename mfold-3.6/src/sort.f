      subroutine sortout(i,j,rep,err)
      include 'rna.inc'
      logical mark
c     The first time in (REP = 1), valid I,J base-pairs are
c     sorted by energy.
      err = 0
      if (rep.eq.1.or.cntrl(7).eq.2) then
         call build_heap
         call heap_sort
c         write(45,*) 'Sorted values'
c         do k = 1,num
c            write(45,*) heapi(k),heapj(k),ene(heapi(k),heapj(k))
c         enddo
         cntr = num
      endif
c     Select the next valid unmarked base-pair
      do while (mark(heapi(cntr),heapj(cntr)))
         if (cntr.eq.1) then
            err = 30
            return
         endif
         cntr = cntr - 1
      enddo
c     The base-pair I,J will be used to create a folding.
      i = heapi(cntr)
      j = heapj(cntr)
 
      return
      end
 
 
c     Add i,j to heapi and heapj if the best energy of a folding containing
c     i,j is no greater than a given percent ( cntrl(8) ) of the minimum
c     folding energy.
c     7/15/96 - M. Zuker modifies the html version. Minimum energy increment 
c     is 1 kcal; maximum is 12 kcal
c     2/7/99 - M. Zuker modifies program. If cntrl(8) < 0, then this is
c     the energy increment in 10ths or 100ths of a kcal/mol
c     Free energies are in 10ths of a kcal/mol or 100ths of a kcal/mol

      subroutine build_heap
      include 'rna.inc'
 
      iprec = ifix(prec)
      if (cntrl(8).lt.0) then
         crit = vmin - cntrl(8)
      else   
         einc = ((abs(vmin) + 5)*cntrl(8))/100
         if (einc.lt.iprec.and.usage.eq.'html') then
            crit = vmin + iprec
         elseif (einc.gt.12*iprec.and.usage.eq.'html') then
            crit = vmin + 120*iprec
         else
            crit = vmin + einc
         endif
      endif   
 
      num = 0
      do diag = 3,2*n-1
         if (diag.le.n+1) then
            i = 1
            j = diag - i
         else
            i = diag - n
            j = diag - i
         endif
         oldene = infinity
         do while (j.gt.i)
            if (ene(i,j).le.crit.and.ene(i,j).ne.oldene) then
               if (num.eq.sortmax) then
                  err = 31
                  call errmsg(err,hstnum(i),hstnum(j))
                  goto 10
               endif
               num =  num + 1
               oldene = ene(i,j)
               heapi(num) = i
               heapj(num) = j
               i = i + 1
               j = j - 1
            else
               i = i + 1
               j = j - 1
            endif
         enddo
      enddo

c      write(45,*) 'crit = ',crit
c      do k = 1,num
c         write(45,*) heapi(k),heapj(k),ene(heapi(k),heapj(k))
c      enddo
      do i = num+1,sortmax+1
        heapi(i) = 0
        heapj(i) = 0
      enddo
 
10    do q = 2,num
         cur = q
         up = cur/2
         do while (up.ge.1.and.ene(heapi(cur),heapj(cur)).lt.ene(heapi(up),heapj(up)))
            call swap(heapi(cur),heapi(up))
            call swap(heapj(cur),heapj(up))
            cur = cur/2
            up = cur/2
         enddo
      enddo
      return
      end
 
 
c     Efficient sort of heap.
      subroutine heap_sort
      include 'rna.inc'
      
      do ir = num-1,1,-1
         call swap(heapi(ir+1),heapi(1))
         call swap(heapj(ir+1),heapj(1))
         up = 1
         c = 2
         do while (c.le.ir)
            if (c.ne.ir) then
               if (ene(heapi(c+1),heapj(c+1)).lt.ene(heapi(c),heapj(c))) then
                  c = c + 1
               endif
            endif
            if (ene(heapi(c),heapj(c)).lt.ene(heapi(up),heapj(up))) then
               call swap(heapi(c),heapi(up))
               call swap(heapj(c),heapj(up))
               up = c
               c = 2 * c
            else
               c = ir+1
            endif
         enddo
      enddo
      return
      end
 
c     ene(k) is the minimum energy of a folding containing the base-pair
c     i,j at heap(k).
      function ene(i,j)
      include 'rna.inc'
      ene = v(i,j) + v(j,i+n)
      return
      end
 
      subroutine swap(i,j)
      k = i
      i = j
      j = k
      return
      end
