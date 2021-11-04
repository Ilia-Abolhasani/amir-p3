c______________________________________________________________________________
c     mfold version 3.6
c     Michael Zuker, Rensselaer Polytechnic Institute
c------------------------------------------------------------------------------
c     - Prediction of RNA/DNA secondary structure by free energy 
c     minimization.

      include 'rna.inc'
      character*30 file_name
      real energy
      logical flag,mark
c     Command line arguments: 
c     1. l for a linear molecular (default) and c for a circular molecule
c     2. text or html
      if (iargc().eq.0) then
         lorc = 'l'
         usage = 'text'
      elseif (iargc().eq.1) then
         call getarg(1,lorc)
         if (lorc.eq.'-v'.or.lorc.eq.'-V') then
            write (6,*) 'nafold from mfold version 3.6'
            call exit(0)
         endif
         usage = 'text'
      else
         call getarg(1,lorc)
         call getarg(2,usage)
      endif
      if (lorc.ne.'l') lorc = 'c'
      if (usage.ne.'text') usage = 'html'

c     Display information
      write(6,*) 'mfold version 3.6 by Michael Zuker'

c     Read file name (not used yet, give CR to continue)
c     file_name = '                              '
c     write(6,*) 'Enter file name: '
c     read(5,2,end=999) file_name
c2    format(a30)

c     Initial setup for run.
      cntrl(3) = 0
5     call device
      if (cntrl(1).ne.2) then
c        Read energy information if this is not a continuation run.
         call enefiles
         call ergread
      else
c        dump out information read in from continuation file
         call cdump
      endif
c     Determine output specifications if this is not a save run.
      if (cntrl(1).ne.1) call outputs
c     Call the menu if this is not a continuation run.
      if (cntrl(1).ne.2) call menu
      mrep = 1
c     CNTRL(7) = 0  -  suboptimal dot plot
c                1  -  N best sorted by energy
c                2  -  best folding for all sequences in a file
10    if (cntrl(7).eq.2) call mseq(mrep)
c     Process sequence before folding.
      do i = 1,mxbits
         marks(i) = 0
         force2(i) = 0
         force2(n+i) = 0
      enddo
      call process
c     Write out processed sequence fragment with auxiliary information
c     underneath. 
      nrows = (n-1)/50 + 1
      do j = 1,nrows
         last = min0(n,50*j)
         lastn = 10*((last + 9)/10)
         write(4,100) (i,i=50*(j-1)+10,lastn,10)
 100     format(5(1x,i10))
         write(4,102) (seq(nsave(1)+i-1),i=1+50*(j-1),last)
         write(4,102) (aux(i),i=1+50*(j-1),last)
 102     format(5(1x,10a1))
      enddo
      if (n*2.gt.fldmax) then
         tt = fldmax/2
c        Fragment is too long. Try again.
         write(6,*) 'STOP: Segment larger than ',tt
         call exit(1)
      endif
      if (cntrl(1).ne.2) then
c        Fill the optimal energy arrays except in a continuation run.
         nofold = 0
         call fill(nofold)
         if (nofold.eq.1) goto 910
      endif
 
      if (cntrl(1).eq.1) then
c        Save the results from FILL in a SAVE run and then stop.
         call putcont
         call exit(0)
      endif
 
      rep = 1
      flag = .true.
      err = 0
      
      do while (flag)
         if (cntrl(7).eq.0) then
c     Interactive dot plot returns iret, jret (new numbering).
            write(6,*) 'Interactive dot plot disabled in this version'
            write(6,*) 'Enter base pair for traceback'
            read(5,*,end=999) iret,jret
            iret = newnum(iret)
            jret = newnum(jret)
            write(6,*) 'V(i,j)= ',v(iret,jret),' V(j,i)= ',v(jret,iret+n)
            write(6,*) 'E(i,j)= ',v(iret,jret) + v(jret,iret+n)
            write(6,*) 'W(i,j)= ',w(iret,jret),' W(j,i)= ',w(jret,iret+n)
            write(6,*) 'W5(i-1)= ',w5(iret-1)
            write(6,*) 'W3(j+1)= ',w3(jret+1)
         else
c     Automatic sort returns IRET,JRET (new numbering).
            call sortout(iret,jret,rep,err)
            if (err.eq.30) then
               flag = .false.
               call errmsg(err,rep-1,0)
               err = 0
            endif
            rep = rep + 1
         endif
c     First traceback yields the best structure on the included fragment
c     from IRET to JRET.
         if (flag) call trace(iret,jret,pen1,err)
         if (err.ne.0) call errmsg(err,iret,jret)
         if (flag) then
            it = iret+n
c     Second traceback yields the best structure on the excluded fragment
c     from IRET to JRET.
            call trace(jret,it,pen2,err)
            if (err.ne.0) then
               call errmsg(err,jret,it)
            else
c     The energy of the best structure containing the base-pair IRET,
c     JRET is the sum of the energies of the optimal foldings on
c     the included and excluded fragments.
               ene = v(iret,jret) + v(jret,iret+n) 
               energy = float(ene - cntrl(3)*(pen1 + pen2)) / prec
c     Count the number of new base pairs not within WINDOW
c     of existing base pairs.
               numbp = 0
               do k = 1,n
                  if(k.lt.basepr(k)) then
                     if(.not.mark(k,basepr(k))) numbp = numbp + 1
                  endif
               enddo
               do k = 1,n
                  if(k.lt.basepr(k)) then
c     Mark "traced-back" base pairs and also base-pairs
c     which are close (within WINDOW = CNTRL(9) ).
                     call smark(k,basepr(k))
                     if(cntrl(9).gt.0) then
                        do k1 = -cntrl(9),cntrl(9)
                           do k2 = -cntrl(9),cntrl(9)
                              if(k+k1.gt.0.and.k+k1.lt.basepr(k)+k2.and.
     1                             basepr(k)+k2.le.n) then
                                 call smark(k+k1,basepr(k)+k2)
                              endif
                           enddo
                        enddo
                     endif
                  endif
               enddo
               if(numbp.le.cntrl(9).and.rep.gt.2) then
                  rep = rep - 1
                  go to 900
               endif
               open(19,file='mfold.log',status='unknown')
               write(19,*) rep-1,','
               close(19)
               if (usage.eq.'text') then
                  if (rep.eq.1) then
                     call system('echo "Structure "|tr -cd "[ -z]" > /dev/tty')
                  endif
                  call system('cat mfold.log|tr -cd ",0-9" > /dev/tty')
               endif
               if (usage.eq.'html'.and.
     .              cntrl(2).ne.2.and.cntrl(2).ne.3.and.cntrl(2).ne.6) then
c     Text output
                  if (rep.eq.2) write(cntrl(4),1011)
 1011             format('<html><head><title>Folding energy details</t',
     .                 'itle><META HTTP-EQUIV="Content-Type" CONTENT="',
     .                 'text/html"; charset=UTF-8"></head>',/,
     .                 '<body bgcolor="#dfc890" text="#000000" link=',
     .                 '"#8e2323" vlink="#2f4f4f" alink="#00ffff"><b>')
                  if (rep.lt.10) then
                     write(cntrl(4),1014) rep-1,rep-1
 1014                format('<a name="STRUCTURE_',i1,'"> Structure ',i1,'</a>')
                  elseif (rep.lt.100) then
                     write(cntrl(4),10141) rep-1,rep-1
10141                format('<a name="STRUCTURE_',i2,'"> Structure ',i2,'</a>')
                  else
                     write(cntrl(4),10142) rep-1,rep-1
10142                format(/'<a name="STRUCTURE_',i3,'"> Structure ',i3,'</a>')
                  endif
                  call linout(1,n,energy,iret,jret,err)
                  write(cntrl(4),1015) rep-1
 1015             format('<hr size=3 width=80% noshade><p><!--endof',i4,' -->')
               else if (usage.ne.'html'.and.
     .                 cntrl(2).ne.2.and.cntrl(2).ne.3.and.cntrl(2).ne.6) then
                  if (cntrl(7).ne.2) then
                     write(cntrl(4),1012) rep-1
 1012                format(/'Structure ',i4,/)
                  else
                     write(cntrl(4),10121) mrep,rep-1
10121                format(/'Sequence  ',i4,' Structure ',i4,/)
                  endif
                  call linout(1,n,energy,iret,jret,err)
               endif
               if (err.ne.0) call errmsg(err,iret,jret)
               if (cntrl(2).ge.3.and.cntrl(2).ne.4.and.usage.eq.'html') then
c     Detailed energy output.
                  if (rep.eq.2) write(22,1013)
 1013             format('<html><head><title>Loop Free-Energy Decompos',
     .                 'ition</title><META HTTP-EQUIV="Content-Type" C',
     .                 'ONTENT="text/html"; charset=UTF-8"></head>',/,
     .                 '<body bgcolor="#dfc890" text="#000000" link=',
     .                 '"#8e2323" vlink="#2f4f4f" alink="#00ffff"><h3>',
     .                 /,'<a name="TOP"><u>Loop Free-Energy Decomposit',
     .                 'ion</u></a></h3><b>')
                  if (rep.le.10) then
                     write(22,1014) rep-1,rep-1
                  elseif (rep.le.100) then
                     write(22,10141) rep-1,rep-1
                  else
                     write(22,10142) rep-1,rep-1
                  endif
                  call erg_det_html(ene,iret,jret,rep-1)
                  write(22,1015) rep-1
               elseif (cntrl(2).ge.3.and.cntrl(2).ne.4.and.usage.eq.'text') 
     .                 then
                  write(22,1012) rep-1
                  call erg_det(ene,iret,jret)
               endif
               
               if (mod(cntrl(2),2).eq.0.or.cntrl(2).eq.7) then
c     ct file output.
                  call ct(energy)
               endif
            endif
         endif
         if (cntrl(7).ge.1.and.rep.gt.cntrl(6)) flag = .false.
 900     continue
      enddo
      if (usage.eq.'text') then
         call system('echo " " > /dev/tty')
      endif
c     
c     Multiple sequence option (cntrl(7) = 2)
c     If sequence number (mrep) is < total number of sequences
c     (cntrl(5)), go get another sequence.
c     
 910  if (cntrl(7).eq.2.and.mrep.lt.cntrl(5)) then
         mrep = mrep + 1
         goto 10
      endif
      if (usage.eq.'html'.and.cntrl(2).ne.2.and.cntrl(2).ne.3.
     .     and.cntrl(2).ne.6) then
c     Finish html text output
         write(cntrl(4),1017)
 1017    format('<address><code>mfold</code> version 3.6<br>M. ',
     .        'Zuker, Rensselaer Polytechnic Institute</address>',
     .        '</body></html>')
      endif
 999  call exit(0)
      end
