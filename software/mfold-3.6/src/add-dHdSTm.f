c     Read .out, .out.html, .det or .det.html file and a corresponding .dh file
c     For each "dG =" line in the first input file, read the dG
c     Read the corresponding dH in the .dh file (stop for discrepancies)
c     Compute dS and Tm based on simple melt.
c     Add dH, dS  and Tm to "dG = " line
c     Output is same as first input, but with extra information on 
c     the "dG =" lines.

      implicit real (a-h,o-z)
      integer amax
      character*4 mode
      character*512 infile1,infile2,outfile
      character*5000 record1,record2
      data amax/5000/

      if (iargc().lt.3) stop 
     .     '2 input files, temperature & mode are required'
      call getarg(1,infile1)
      call getarg(2,infile2)
      call getarg(3,record1)
      if (iargc().eq.3) then
         mode = 'text'
      else
         call getarg(4,mode)
      endif
      read(record1,*) t
      t = t + 273.15
      labstart = index(infile2,'.dh')
      outfile = infile2
      outfile(labstart:labstart+6) = '.dHdSTm'
      open(2,file=infile1,status='old',err=91)
      open(3,file=infile2,status='old',err=92)
      open(4,file=outfile,status='unknown')
      
      read(2,100,end=93) record1
 100  format(a5000)
      do while (1.eq.1)
         do while (index(record1,' dG =').eq.0)
            call dump(amax,record1)
            read(2,100,end=99) record1
         enddo
         read(record1(index(record1,' = ')+3:amax),*) dg
c     Now look for dh in second input file.
         read(3,100,end=94) record2
         do while (index(record2,'Computed energy = ').eq.0)
            read(3,100,end=94) record2
         enddo
         read(record2(index(record2,' = ')+3:amax),*) dh

c     Now compute. Round dg to 2 decimals.
         idg = nint(100.0*dg)
         dg = float(idg)/100.0
         ds = (dh - dg)/t
         tm = t*dh/(dh - dg)
         istart = index(record1,' dG =') + 5
         if (mode.eq.'html') then
            write(record1(istart:amax),105) dg,dh,1000.0*ds,tm-273.15,
     .           char(226),char(132),char(131)
 105        format(1x,f9.2,' &nbsp; dH = ',f9.2,' &nbsp; dS = ',f9.2,
     .           ' &nbsp; T<sub>m</sub> = ',f6.1,' ',3a1)
         else
            write(record1(istart:amax),107) dg,dh,1000.0*ds,tm-273.15,
     .           char(226),char(132),char(131)
 107        format(1x,f9.2,'  dH = ',f9.2,'  dS = ',f9.2,'  Tm = ',f6.1,
     .           ' ',3a1)
         endif
         call dump(amax,record1)
         record1 = '                                         '
      enddo

 91   stop 'Error opening first input file'
 92   stop 'Error opening second input file'
 93   stop 'First input file is empty'
 94   stop 'Missing information in dh file'

 99   call exit(0)
      end
      subroutine dump(amax,record)
      implicit real (a-h,o-z)
      integer amax
      character*1 record(amax)
      m = amax
      do while (record(m).eq.' '.and.m.gt.1)
         m = m - 1
      enddo
      write(4,100) (record(l),l=1,m)
 100  format(5000a1)
      return
      end
