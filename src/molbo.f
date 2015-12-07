c-----------------------------------------------------------------------
c--- MOLBO: a utility to generate the *.47 file from MOLPRO output.
c--- Written in FORTRAN 77 and a bit of FORTRAN 90.
c---
c--- Tested for MOLPRO 2006~2010 and NBO 3, 5, and 6.
c---
c--- E-mail: qcband@gmail.com
c-----------------------------------------------------------------------

      program MOLBO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(nmt=8)
      dimension mlst(nmt)

      character*10 dt
      character*5 ver
      character*27 f31,f47
      character*32 fmb

c-----------------------------------------------------------------------
c--- 0. initial
c-----------------------------------------------------------------------
      ver="2.1.5"
      dt="07/04/2014"
      call head(ver,dt)

c--- i/o ports
      imop=45
      i47f=55
      i31f=56
      imbf=57

c--- matrix list
c---  the following matrices are necessary for NBO/MBO analysis
c      mlst(1) = 1    ! overlap
c      mlst(2) = 1    ! density
c---  the following matrices are optional
c      mlst(3) = 0/1  ! fock
c      mlst(4)        ! contract; will be redefined in plotnbo
c              = -2,...+2
c      mlst(5) = 0/1  ! MO
c      mlst(6) = 0/1  ! kinetic
c      mlst(7) = 0/1  ! potential
c      mlst(8) = 0/1  ! dipole(x,y,z)
c---  special
c---   mlst(4) = +/-1 ! spherical
c---           = +/-2 ! cartesian
c---           < 0    ! for NBO3 + plot + G-function

c-----------------------------------------------------------------------
c--- 1. define MOLPRO output, *.31, *.47, and mbo files
c-----------------------------------------------------------------------
      call molout(f31,f47,fmb)

c-----------------------------------------------------------------------
c--- 2. check
c-----------------------------------------------------------------------
c--- 2.1. search dummy atoms
      call rddummy(ist)
      if(ist.eq.1)then
        write(*,"(' Wrong! Dummy atom is found.')")
        goto 9910
      end if

c--- 2.2. check density fitting
      call rddf(idf)
      if(idf.ne.0)write(*,"(' Density fitting calculation.',/)")

c--- 2.3. check DKH
      call rddkh(idk)
      if(idk.ne.0)write(*,"(' DKH relativistic calculation.',/)")

c-----------------------------------------------------------------------
c--- 3. read coordinates
c-----------------------------------------------------------------------
      call rdcoord(ist,idf)
      if(ist.eq.0)then
        write(*,"(' Wrong! Can not find the coordinates!')")
        goto 9910
      else if(ist.eq.-1)then
        write(*,"(' Wrong! Too many atoms! Please inctrease MAXNATM.')")
        goto 9910
      else if((ist.ne.1 .and. idf.eq.0) .or.
     *        (ist.gt.2 .and. idf.eq.1) )then
        write(*,"(' Wrong! ',i3,' groups of coordinates are found!',/,
     *  ' Only single-point calculation is supported.')")ist
        goto 9910
      end if

c-----------------------------------------------------------------------
c--- 4. read basis sets
c-----------------------------------------------------------------------
c--- 4.1. number of contractions
      call rdnbas(nbas,ist,idf)
      if(nbas.eq.0)then
        write(*,"(' Wrong! Can not find the basis set information!')")
        goto 9910
      else if(nbas.lt.0)then
        write(*,"(' Wrong! ',i3,' groups of basis set are found!',/,
     *  ' The basis set can be defined only once.')")abs(nbas)
        goto 9910
      end if
      if(ist.eq.1)then
        write(*,"(
     *  ' Wrong! The number of basis functions is too large!')")
        goto 9910
      end if

c--- 4.2. plot orbitals or not
      call plotnbo(mlst(4))
      if(mlst(4).ne.0)then
        call carsph(mlst(4))
c---  mlst(4) = 0 (unknown; no $CONTRACT part), 1 (spherical), 2 (cartesian)
        if(mlst(4).eq.0)then
          write(*,"(/,' Unknown basis function type.',/,
     *    ' Orbitals can not be plotted.')")
        else if(mlst(4).eq.1)then
          write(*,"(/,' Spherical basis functions are used.')")
        else if(mlst(4).eq.2)then
          write(*,"(/,' Cartesian basis functions are used.')")
        end if
      end if

c--- 4.3. read the basis function types
      call rdbstyp(ist,nbas,ig)
      if(ist.eq.-1)then
        write(*,"(//,' Wrong! Can not find the basis function types!',/,
     *  ' Please use the commond gprint,basis=10 in the calculation.')")
        goto 9910
      else if(ist.eq.1)then
        write(*,"(//,' Wrong! Higher functions than G are found!',/,
     *  ' They are not supported by NBO program.')")
        goto 9910
      else if(ist.eq.2)then
        write(*,"(//,
     *  ' Wrong! Symmetry equivalent atoms can not be treated!',//,
     *  ' Please decrease the symmetry in the MOLPRO calculation',/,
     *  ' until all the atoms are symmetry unique.')")
        goto 9910
      end if

c--- 4.4. reorder basis functions according to atoms and L-quantum numbers
c---    this is necessary for plot
c--- DEBUG
c      if(mlst(4).ne.0) call robas(nbas)
c---
      call robas(nbas)

c--- 4.5. read basis functions and contractions
c---    this is necessary for "plot"
      if(mlst(4).ne.0)then
        call rdfun(mlst(4),nbas)
        if(mlst(4).eq.0)then
          write(*,"(/,' No basis functions are read,'
     *' so orbitals can not be plotted.',//)")
        else if(mlst(4).eq.-100)then
          write(*,"(/,' Wrong! sum(NCOMP) .ne. NBAS,'
     *' so orbitals can not be plotted.',//)")
          mlst(4)=0
        else if(mlst(4).eq.-200)then
          write(*,"(/,' Wrong! sum(NPRIM) .ne. NEXP,'
     *' so orbitals can not be plotted.',//)")
          mlst(4)=0
        end if
      end if
      write(*,"(1x,45('-'),//)")

c-----------------------------------------------------------------------
c--- 5. the number of symmetry blocks
c-----------------------------------------------------------------------
      call rdsymm(ist)
      if(ist.eq.0)then
        write(*,"(' Wrong! Can not find the point group!')")
        goto 9910
      end if
c--- # b.s. in each block
      call rdblk(nbas,ist)
      if(ist.eq.0)then
        write(*,"(' Wrong! NBAS .ne. NETOT!')")
        write(*,"(' NBAS=',i4,' NETOT=',i4)")nbas,netot
        goto 9910
      end if

c-----------------------------------------------------------------------
c--- 6. overlap matrix
c--- mlst(1) =-1  wrong format
c---         = 0  no overlap matrix
c---         = 1  normal
c---         = 2  use lower symmetry
c---         = 3  unknown reason for C1 symmetry
c-----------------------------------------------------------------------
      call rdovlp(mlst(1))
      if(mlst(1).eq.0)then
        write(*,"(' *** Error',/,
     *  ' No overlap matrix is found, which is necessary.',//)")
        goto 9910
      else if(mlst(1).eq.-1)then
        write(*,"(' Error when reading the overlap matrix!')")
        goto 9910
      else if(mlst(1).ne.1)then
        write(*,"(
     *' Wrong diagonal elements of the overlap matrix since the matrix',
     */,' has been block-diagonalized by MOLPRO.'
     *  )")
        if(mlst(1).eq.2)then
          write(*,"(/,
     *' In order to avoid this problem, please decrease the symmetry in'
     *,/,' MOLPRO until all the atoms being symmetry unique.')")
        else
          write(*,"(/,
     *' Unknown problem for the C1 point group.')")
        end if
        goto 9910
      end if

c-----------------------------------------------------------------------
c--- 7. density matrix
c-----------------------------------------------------------------------
      call rddens(mlst(2))
      if(mlst(2).eq.0)then
        write(*,"(' *** Error',/,
     *  ' No density matrix is found, which is necessary.',//)")
        goto 9910
      else if(mlst(2).eq.-1)then
        write(*,"(' Error when reading the density matrix!')")
        goto 9910
      end if

c-----------------------------------------------------------------------
c--- 8. Fock matrix; this matrix is optional
c-----------------------------------------------------------------------
      call rdfock(ist)
      if(ist.eq.0)then
        mlst(3)=0
        write(*,"(' *** Warning',/,
     *  ' No Fock matrix is found, which is optional.',//)")
      else if(ist.eq.-1)then
        write(*,"(' Error when reading the Fock matrix!')")
        goto 9910
      else
        mlst(3)=1
      end if

c-----------------------------------------------------------------------
c--- 9. MO matrix; this matrix is optional
c-----------------------------------------------------------------------
      call rdmo(ist)
      if(ist.eq.0)call rdmo2(ist)
      if(ist.eq.0)then
        mlst(5)=0
        write(*,"(' *** Warning',/,
     *  ' No MO matrix is found, which is optional.',//)")
      else if(ist.eq.-1)then
        write(*,"(' Error when reading the MO matrix!')")
        goto 9910
      else
        mlst(5)=1
      end if

c-----------------------------------------------------------------------
c--- 10. kinetic matrix; this matrix is optional
c-----------------------------------------------------------------------
      if(idk.eq.0)then
        call rdkine(ist)
        if(ist.eq.0)then
          mlst(6)=0
          write(*,"(' *** Warning',/,
     *    ' No kinetic matrix is found, which is optional.',//)")
        else if(ist.eq.-1)then
          write(*,"(' Error when reading the kinetic matrix!')")
          goto 9910
        else
          mlst(6)=1
        end if
      else
        write(*,
     *   "(' Non-relativistic kinetic matrix is omitted by DKH.',/)")
      end if

c-----------------------------------------------------------------------
c--- 11. potential matrix; this matrix is optional
c-----------------------------------------------------------------------
      if(idk.eq.0)then
        call rdpote(ist)
        if(ist.eq.0)then
          mlst(7)=0
          write(*,"(' *** Warning',/,
     *    ' No potential matrix is found, which is optional.',//)")
        else if(ist.eq.-1)then
          write(*,"(' Error when reading the potential matrix!')")
          goto 9910
        else
          mlst(7)=1
        end if
      else
        write(*,
     *   "(' Non-relativistic potential matrix is omitted by DKH.',/)")
      end if

c-----------------------------------------------------------------------
c--- 12. dipole matrix (x,y,z); this matrix is optional
c-----------------------------------------------------------------------
      call rddip(ist)
      mlst(8)=0
      if(ist.eq.0)then
        write(*,"(/,
     *  ' Symmetry is used, so the $DIPOLE datalist is OFF.',//)")
      else if(ist.eq.-1)then
        write(*,"(' *** Warning',/,
     *  ' No dipole matrix (x) is found, which is optional.',//)")
      else if(ist.eq.-2)then
        write(*,"(' Error when reading the dipole matrix (x)!')")
        goto 9910
      else if(ist.eq.-3)then
        write(*,"(' *** Warning',/,
     *  ' No dipole matrix (y) is found, which is optional.',//)")
      else if(ist.eq.-4)then
        write(*,"(' Error when reading the dipole matrix (y)!')")
        goto 9910
      else if(ist.eq.-5)then
        write(*,"(' *** Warning',/,
     *  ' No dipole matrix (z) is found, which is optional.',//)")
      else if(ist.eq.-6)then
        write(*,"(' Error when reading the dipole matrix (z)!')")
        goto 9910
      else if(ist.eq.1)then
        mlst(8)=1
      end if

c-----------------------------------------------------------------------
c--- 13. NBO version. inbo: 0=new; 1=old
c-----------------------------------------------------------------------
      if(mlst(4).ne.0 .and. ig.ne.0)then
        write(*,"(' G-functions are found for plotting.')")
        call nbover(inbo)
c       G-functions are not supported by NBO 3
        if(inbo.eq.1) call file31(mlst(4))
      end if

c-----------------------------------------------------------------------
c--- 14. write *.47, *.31 for NBO
c-----------------------------------------------------------------------
      OPEN(i47f,FILE=f47)
      rewind(i47f)
      if(mlst(4).lt.0)then
        OPEN(i31f,FILE=f31)
        rewind(i31f)
      end if

      call wrnbokey(nbas,mlst(4))
      call wrcoord(ig,mlst(4))
      call wrbas(nbas)
      call wrctr(mlst,ig)
      call wrmat(mlst)

c-----------------------------------------------------------------------
c--- 15. Mayer's Bond Orders
c-----------------------------------------------------------------------
      OPEN(imbf,FILE=fmb)
      rewind(imbf)
      call mayer

c-----------------------------------------------------------------------
c--- 16. the final step
c-----------------------------------------------------------------------
      write(*,"(///,1x,45('-'),/,2x,
     *'Following files are generated successfully!',/,1x,45('-'))")
      if(mlst(4).lt.0)then
        write(*,"(2x,a)")trim(f31)
        close (i31f)
      end if
      write(*,"(2x,a,/,2x,a)")trim(f47),trim(fmb)
      close (i47f)
      close (imbf)
      call datalist

9910  continue
      close (imop)

      call estop

      end

c-----------------------------------------------------------------------
c--- version of NBO
c--- 0=new; 1=old
c-----------------------------------------------------------------------
      subroutine nbover(inbo)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*1 ver

100   write(*,"(/,
     *' Which version of NBO program do you use?',/,
     *'  1) NBO 4 or higher (default)',/,
     *'  2) NBO 3',/,
     *' > ',$)")
      read(*,"(a1)")ver
      lenth=len_trim(ver)
      if(lenth.eq.0.or.ver.eq.'1')then
        inbo=0
      else if(ver.eq.'2')then
        inbo=1
      else
        write(*,1000)ver
        goto 100
      end if

      return
1000  format(' Unknown "',a1,'". Try again.')
      end

c-----------------------------------------------------------------------
c--- Mayer's Bond Orders
c---   MBO(I,J) = sum_iI (sum_iJ ( (PxS)_iJ,iI * (PxS)_iI,iJ ) )
c--- Since alpha and beta densities are merged in the output, the
c--- calculated MBOs are equivalent to the generalized Wiberg bond order
c--- indices (for AO basis instead of NBO basis).
c--- Reference
c---   S. I. Gorelsky, AOMix manual.
c--- Generalized Wiberg bond order indices
c---   Mayer, I. Chem. Phys. Lett. 1983, 97, 270.
c--- MBO
c---   Mayer, I. Chem. Phys. Lett. 1984, 110, 440.
c---   Mayer, I. Theor. Chim. Acta 1985, 67, 315.
c---   Mayer, I. Int. J. Quantum Chem. 1986, 29, 73.
c---   Mayer, I. Int. J. Quantum Chem. 1986, 29, 477.
c-----------------------------------------------------------------------
      subroutine mayer
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000,maxnatm=100,max_za=103)
      common/COORD/xyz(3,maxnatm),iz(2,maxnatm),natm
      common/BASIS/icent(maxbas),label(maxbas)
      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      dimension ps(maxbas,maxbas),boao(maxbas,maxbas),
     *bo(maxnatm,maxnatm)
      character*2 ATOMLIB(max_za)
      data (atomlib(i),i=1,50) /
     1'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     2'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     3'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     4'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     5'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn'/
      data (atomlib(i),i=51,100) /
     1'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     2'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     3'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     4'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     5'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/
      data (atomlib(i),i=101,103) /
     1'Md','No','Lr'/

      imbf=57

      do j=1,netot
        do i=1,netot
          ps(i,j)=0.d0
          do k=1,netot
            ps(i,j)=ps(i,j)+dens(i,k)*ovlp(k,j)
          end do
        end do
      end do

c      do i=1,netot
c        write(9,'(100f10.6)')(ps(j,i), j=1,netot)
c      end do

      do j=1,netot
        do i=1,netot
          boao(i,j)=ps(i,j)*ps(j,i)
        end do
      end do

c      do i=1,netot
c        write(9,'(100f10.6)')(boao(j,i), j=1,netot)
c      end do

c--- DENS and OVLP have been reordered in the subroutine WRMAT
      bo=0.d0
      do i=1,netot
        do j=i+1,netot
          ii=icent(i)
          jj=icent(j)
          if(ii.ne.jj)then
            bo(ii,jj)=bo(ii,jj)+boao(i,j)
c            bo(jj,ii)=bo(ii,jj)    ! not used
            bo(ii,ii)=bo(ii,ii)+boao(i,j)
            bo(jj,jj)=bo(jj,jj)+boao(i,j)
          end if
        end do
      end do

      write(imbf,"(' Mayer bond order matrix:',/)")
      write(imbf,"(4x,'Atom',100(2x,i3,'-',a2))")
     *(i,atomlib(iz(1,i)),i=1,natm)
      do i=1,natm
        write(imbf,"(2x,i3,'-',a2,100f8.4)")i,atomlib(iz(1,i)),
     *    (bo(j,i), j=1,i)
      end do

      write(imbf,"(/,
     *' E(i,i) = Total MBO of atom i, and',/,
     *' E(i,j) = MBO between atoms i and j.')")

      return
      end

c-----------------------------------------------------------------------
c--- generate the .31 file
c-----------------------------------------------------------------------
      subroutine file31(ifplot)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*1 yn

      write(*,"(/,
     *' G-functions can not be used in the .31 file by NBO 3.',/,
     *' Do you want to generate the .31 file by MOLBO? ([Y] or N)')")
100   write(*,"(' > ',$)")
      read(*,"(a1)")yn
c--- default
      if(len_trim(yn).eq.0)yn='Y'

      if(yn.eq.'Y'.or.yn.eq.'y')then
        ifplot=-ifplot
      else if(yn.eq.'N'.or.yn.eq.'n')then
        ifplot=0
        write(*,"(/,' The PLOT keyword is OFF.')")
      else
        write(*,"(/,' Unknown selection. Please try again.',/)")
        goto 100
      end if

      return
      end

c-----------------------------------------------------------------------
c--- write $CONTRACT
c-----------------------------------------------------------------------
      subroutine wrctr(mlst,ig)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      common/BASIS/icent(maxbas),label(maxbas)
      parameter(maxexp=5000,maxgroup=800)
      common/CONTRACT/NSHELL,NEXP,NCOMP(maxgroup),NPRIM(maxgroup),   ! for NBO
     *NPTR(maxgroup),EX(maxexp),CS(maxexp),CP(maxexp),CD(maxexp),
     *CF(maxexp),CG(maxexp)
      dimension mlst(*)
      character*10 tag(10)
      data tag/
     *'   NCOMP =','   NPRIM =','    NPTR =','          ','     EXP =',
     *'      CS =','      CP =','      CD =','      CF =','      CG ='/

      i47f=55
      i31f=56

      if(mlst(4).eq.0)goto 9999

      if(mlst(4).gt.0)goto 5000
c--- .31
      idx=0
      do i=1,NSHELL
        write(i31f,"(1x,4i6)")icent(idx+1),NCOMP(i),NPTR(i),NPRIM(i)
        write(i31f,"(1x,10i6)")label(idx+1:idx+NCOMP(i))
        idx=idx+NCOMP(i)
      end do
      write(i31f,"(1x,75('-'))")

      write(i31f,"(2x,4e18.9)")(EX(i),i=1,NEXP)
      write(i31f,"()")
      write(i31f,"(2x,4e18.9)")(CS(i),i=1,NEXP)
      write(i31f,"()")
      write(i31f,"(2x,4e18.9)")(CP(i),i=1,NEXP)
      write(i31f,"()")
      write(i31f,"(2x,4e18.9)")(CD(i),i=1,NEXP)
      write(i31f,"()")
      write(i31f,"(2x,4e18.9)")(CF(i),i=1,NEXP)
      write(i31f,"()")
      write(i31f,"(2x,4e18.9)")(CG(i),i=1,NEXP)
      goto 9999

c--- .47
5000  continue
      write(i47f,"(' $CONTRACT')")
      write(i47f,"('  NSHELL =',i4)")NSHELL
      write(i47f,"('    NEXP =',i4)")NEXP

      nstep=17
      call lines(NSHELL,nstep,nline,last)
c NCOMP
      call iwrite(0,tag(1),tag(4),i47f,nline,last,nstep,NCOMP)
c NPRIM
      call iwrite(0,tag(2),tag(4),i47f,nline,last,nstep,NPRIM)
c NPTR
      if(NPTR((nline-1)*nstep+last).lt.1000)then
        call iwrite(0,tag(3),tag(4),i47f,nline,last,nstep,NPTR)
      else
        call iwrite(1,tag(3),tag(4),i47f,nline,last,nstep,NPTR)
      end if

      nstep=4
      call lines(NEXP,nstep,nline,last)
c EXP
      call fwrite(tag(5),tag(4),i47f,nline,last,nstep,EX)
c CS
      call fwrite(tag(6),tag(4),i47f,nline,last,nstep,CS)
c CP
      call fwrite(tag(7),tag(4),i47f,nline,last,nstep,CP)
c CD
      call fwrite(tag(8),tag(4),i47f,nline,last,nstep,CD)
c CF
      call fwrite(tag(9),tag(4),i47f,nline,last,nstep,CF)
c CG
      if(ig.ne.0)call fwrite(tag(10),tag(4),i47f,nline,last,nstep,CG)

      write(i47f,"(' $END')")

9999  continue
      return

      end

c-----------------------------------------------------------------------
c--- write with the format a, n*i
c-----------------------------------------------------------------------
      subroutine iwrite(Mode,tag1,tag2,io,nline,last,nstep,ivec)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*10 tag1,tag2,ac
      dimension ivec(*)

      idx=0
      if(Mode.eq.0)then
        do i=1,nline
          if(i.eq.1)then
            ac=tag1
          else
            ac=tag2
          end if
          if(i.lt.nline)then
            write(io,1110)ac,(ivec(j),j=idx+1,idx+nstep)
            idx=idx+nstep
          else
            write(io,1110)ac,(ivec(j),j=idx+1,idx+last)
          end if
        end do
      else
        do i=1,nline
          if(i.eq.1)then
            ac=tag1
          else
            ac=tag2
          end if
          if(i.lt.nline)then
            write(io,1111)ac,(ivec(j),j=idx+1,idx+nstep)
            idx=idx+nstep
          else
            write(io,1111)ac,(ivec(j),j=idx+1,idx+last)
          end if
        end do
      end if

      return

1110  format(a10,17i4)
1111  format(a10,17i5)
      end

c-----------------------------------------------------------------------
c--- write with the format a, n*f
c-----------------------------------------------------------------------
      subroutine fwrite(tag1,tag2,io,nline,last,nstep,fvec)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*10 tag1,tag2,ac
      dimension fvec(*)

      idx=0
      do i=1,nline
        if(i.eq.1)then
          ac=tag1
        else
          ac=tag2
        end if
        if(i.lt.nline)then
          write(io,1110)ac,(fvec(j),j=idx+1,idx+nstep)
          idx=idx+nstep
        else
          write(io,1110)ac,(fvec(j),j=idx+1,idx+last)
        end if
      end do

      return

1110  format(a10,4e16.7)
      end

c-----------------------------------------------------------------------
c--- cartesian vs spherical basis functions
c---  itp = 0 (unknown), 1 (spherical), 2 (cartesian)
c-----------------------------------------------------------------------
      subroutine carsph(itp)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*32 tmp,tags1,tagc1
      character*16 tags2,tagc2

      imop=45
      itp=0
      tags1=' Using spherical harmonics'
      tagc1=' Using cartesian basis functions'
      tags2=' Sphericals:   T'
      tagc2=' Sphericals:   F'

      rewind(imop)
10    read(imop,"(a32)",end=300)tmp
      if(index(tmp,tags1)+index(tmp,tags2).ne.0)then
        itp=1
        goto 300
      else if(index(tmp,tagc1)+index(tmp,tagc2).ne.0)then
        itp=2
        goto 300
      end if
      goto 10

300   return
      end

c-----------------------------------------------------------------------
c--- read basis functions
c-----------------------------------------------------------------------
      subroutine rdfun(ibstyp,nbas)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxexp=5000,maxgroup=800,maxctra=25,maxprim=40)
      common/CONTRACT/NSHELL,NEXP,NCOMP(maxgroup),NPRIM(maxgroup),   ! for NBO
     *NPTR(maxgroup),EX(maxexp),CS(maxexp),CP(maxexp),CD(maxexp),
     *CF(maxexp),CG(maxexp)
      character*53 tmp,tagend
      character*27 tagstt
      character*17 tagini
      character*5 cngroup
      dimension nqm(maxgroup),lq(maxgroup),expo(maxprim),
     *cf0(maxprim,maxctra),
     *NPRIMS(maxgroup),NCONTS(maxgroup),NFUNC(maxgroup)              ! for MOLPRO

      if(ibstyp.eq.0)goto 9999

      imop=45
      ibstmp=50
      tagini='PROGRAM * SEWARD '
      tagstt=' Number of groups:         '
      tagend=' DATASETS  * FILE   NREC   LENGTH (MB)   RECORD NAMES'
c--- The first TAGEND is after the integration.
c---
c--- Since MOLPRO 2006, there is another TAGEND in DF-HF/DF-DFT for density fitting.
c--- Since MOLPRO 2008, there is another TAGEND in HF/DFT in which smaller BS is used
c--- for ATOMIC-GUESS.
c--- Both of them must be skipped.
c--- For restart calculation from file-2, we should find TAGINI to skip the first TAGEND.

      rewind(imop)
1     read(imop,"(a53)",end=5)tmp
      if(index(tmp,tagini).ne.0) goto 10
      goto 1
5     ibstyp=0
      goto 9999

10    read(imop,"(a53)",end=50)tmp
      if(index(tmp,tagstt).ne.0) goto 100
      if(index(tmp,tagend).ne.0) goto 50
      goto 10
50    ibstyp=0
      goto 9999

100   continue
c read # Number of groups:
      write(cngroup,"(a5)")tmp(33:37)
      read(cngroup,"(i5)")ngroup
      if(ngroup.le.0)then
        ibstyp=0
        goto 9999
      else if(ngroup.gt.maxgroup)then
        write(*,"(' ### Too many groups in the basis set.')")
        write(*,"(2x,'NGROUP=',i8,' MAXGROUP=',i8)")ngroup,maxgroup
        ibstyp=0
        goto 9999
      end if
c      write(*,*)'ngroup=',ngroup

      read(imop,"(////)")
      NSHELL=0
      do i=1,ngroup
        read(imop,"(5x,3i8)")NPRIMS(i),NCONTS(i),NFUNC(i)
        if(NCONTS(i).gt.maxctra)then
          write(*,"(/,' ### Too many contractions.')")
          write(*,"(2x,'NCONTS=',i8,' MAXCTRA=',i8,/)")NCONTS(i),maxctra
          ibstyp=0
          goto 9999
        end if
        if(NPRIMS(i).gt.maxprim)then
          write(*,"(/,' ### Too many primitive basis functions.')")
          write(*,"(2x,'NPRIMS=',i8,' MAXPRIM=',i8,/)")NPRIMS(i),maxprim
          ibstyp=0
          goto 9999
        end if
        NSHELL=NSHELL+NCONTS(i)
        nqm(i)=NFUNC(i)/NCONTS(i)
        lq(i)=lquant(ibstyp,nqm(i))
        if(lq(i).lt.0)then
          write(*,"(/,' Wrong LQ! Skip basis functions.',/)")
          ibstyp=0
          goto 9999
        end if
c        write(*,"(5x,3i8)")NPRIMS(i),NCONTS(i),NFUNC(i)
c        write(*,"(5x,2i8)")nqm(i),lq(i)
      end do
c      write(*,"(i8)")NSHELL

      open(ibstmp,file='bs123456789.tmp',status='new')
      rewind(ibstmp)

      read(imop,"(/)")
      NEXP=0
      idx=0
      do i=1,ngroup
        read(imop,*)
c--- The printed basis set has been re-normalized by MOLPRO, so the
c--- re-normalization is not necessary.
        do j=1,NPRIMS(i)
          read(imop,*)expo(j),(cf0(j,k),k=1,NCONTS(i))
        end do
        do j=1,NCONTS(i)
          idx=idx+1
          npnew=NPRIMS(i)
          do k=1,NPRIMS(i)
            if(abs(cf0(k,j)).lt.1.d-10)npnew=npnew-1
          end do
          NCOMP(idx)=nqm(i)
          NPRIM(idx)=npnew
          NEXP=NEXP+npnew
          write(ibstmp,"(2i4)")lq(i),npnew
          do k=1,NPRIMS(i)
            if(abs(cf0(k,j)).gt.1.d-10)
     *      write(ibstmp,"(2f20.10)")expo(k),cf0(k,j)
          end do
        end do
      end do

      NPTR(1)=1
      do i=2,NSHELL
        NPTR(i)=NPTR(i-1)+NPRIM(i-1)
      end do

c--- read basis functions
      rewind(ibstmp)
      idx=0
      do i=1,NSHELL
        read(ibstmp,"(2i4)")lqnew,npnew
        do j=1,npnew
          idx=idx+1
          read(ibstmp,"(2f20.10)")EX(idx),cfnew
c--- cfnew = cfnew * (normalization factor)
          call fnorm(lqnew,EX(idx),cfnew)
c--- copy to cs/cp/cd/cf/cg
          select case(lqnew)
            case(0)
              CS(idx)=cfnew
            case(1)
              CP(idx)=cfnew
            case(2)
              CD(idx)=cfnew
            case(3)
              CF(idx)=cfnew
            case(4)
              CG(idx)=cfnew
          end select
        end do
      end do

c--- check
c--- sum(NCOMP)==nbas
      n=0
      do i=1,NSHELL
        n=n+NCOMP(i)
      end do
      if(n.ne.nbas)then
        write(*,"(' sum(NCOMP)=',i4,' NBAS=',i4)")n,nbas
        ibstyp=-100
        goto 9990
      end if
c--- sum(NPRIM)==idx==NEXP
      n=0
      do i=1,NSHELL
        n=n+NPRIM(i)
      end do
      if(n.ne.NEXP .or. idx.ne.NEXP)then
        write(*,"(' sum(NPRIM)=',i4,' NEXP=',i4,' IDX=',i4)")n,NEXP,idx
        ibstyp=-200
        goto 9990
      end if

9990  continue
      close(ibstmp,status='delete')

9999  return
      end

c-----------------------------------------------------------------------
c--- calculate the normalization factor for GTO(L)
c-----------------------------------------------------------------------
      subroutine fnorm(lq,ex,cf)
      implicit double precision (a-h,o-z)

      pi=acos(-1.d0)
      pi3=pi**3.d0

c---  unnormalize primitives
c---  Normal^4 = 2^n1 * a^n2 / (pi^3 * nf)
c---  where n1=3+4*L; n2=3+2*L, nf=[(2L-1)!!]^2
      call power(lq,n1,n2,nf)
      f = (2.d0**dble(n1)) * (ex**dble(n2)) / (pi3 * dble(nf))
      cf=cf*sqrt(sqrt(f))

      return
      end

c-----------------------------------------------------------------------
c--- get power(n1,n2,nf) for GTO(L) normalization
c--- n1=3+4*L; n2=3+2*L, nf=[(2L-1)!!]^2
c-----------------------------------------------------------------------
      subroutine power(lq,n1,n2,nf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      n1=0
      n2=0
      nf=0
      select case(lq)
        case(0)  ! s
          n1=3
          n2=3
          nf=1
        case(1)  ! p
          n1=7
          n2=5
          nf=1
        case(2)  ! d
          n1=11
          n2=7
          nf=9
        case(3)  ! f
          n1=15
          n2=9
          nf=225
        case(4)  ! g
          n1=19
          n2=11
          nf=11025
      end select

      return
      end

c-----------------------------------------------------------------------
c--- L-quantum number
c-----------------------------------------------------------------------
      function lquant(ibs,nq)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      l=-1

      if(ibs.eq.0)goto 999

      if(ibs.eq.2)goto 100

c--- spherical
      select case(nq)
        case(1)
          l=0
        case(3)
          l=1
        case(5)
          l=2
        case(7)
          l=3
        case(9)
          l=4
      end select
      goto 999

c--- cartesian
100   continue
      select case(nq)
        case(1)
          l=0
        case(3)
          l=1
        case(6)
          l=2
        case(10)
          l=3
        case(15)
          l=4
      end select

999   lquant=l

      return
      end

c-----------------------------------------------------------------------
c--- reorder basis functions according to atoms and L-quantum numbers
c--- NOTE: The new order will be wrong if there are symmetry equivalent
c--- atoms. In this case, lower symmetry should be used.
c--- This has been checked in subroutine rdbstyp.
c-----------------------------------------------------------------------
      subroutine robas(nbas)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100,maxbas=2000)
      common/COORD/xyz(3,maxnatm),iz(2,maxnatm),natm
      common/BASIS/icent(maxbas),label(maxbas)
      common/RELATION/lmap(maxbas)
      dimension nml(5),itmp(maxbas)
      data nml/1,3,6,10,15/

      itmp=0
      lmap=0

      idx=0
      if1=0
      if2=0

      do i=1,natm
20      if(if1.ne.0)then
          if1=0
          cycle
        else
          if1=1
        end if
        do lq=1,5            ! S,P,D,F,G
40        if(if2.ne.0)then
            if2=0
            cycle
          else
            if2=1
          end if
c---  ic1 for cartesian bs, ic2 for spherical bs
          ic0=100*(lq-1)
          ic1=ic0
          if(lq.gt.2)then    ! D,F,G
            ic2=ic1+50
          else               ! S,P
            ic2=ic1
          end if
          do lm=1,nml(lq)
            ic1=ic1+1
            ic2=ic2+1
            do j=1,nbas
              if(icent(j).eq.i)then
                if(itmp(j).eq.1)cycle
                if(label(j).eq.ic1.or.label(j).eq.ic2)then
                  idx=idx+1
                  lmap(idx)=j
                  itmp(j)=1
                  if1=0
                  if2=0
                  exit
                end if
              end if
            end do
          end do
          goto 40
        end do
        goto 20
      end do

c--- DEBUG
c      write(8,"(10i4)")icent(1:nbas)
c      write(8,"(10i4)")label(1:nbas)
c      write(8,"(10i4)")lmap(1:nbas)
c      write(8,"(10i4)")icent(lmap(1:nbas))
c      write(8,"(10i4)")label(lmap(1:nbas))
c---
      do i=1,nbas
        itmp(i)=icent(lmap(i))
      end do
      icent(1:nbas)=itmp(1:nbas)
      do i=1,nbas
        itmp(i)=label(lmap(i))
      end do
      label(1:nbas)=itmp(1:nbas)

      return
      end

c-----------------------------------------------------------------------
c--- check DKH
c-----------------------------------------------------------------------
      subroutine rddkh(idkh)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*34 tmp,tag

      imop=45
      tag=' Computing Douglas-Kroll integrals'

      idkh=0
      rewind(imop)
10    read(imop,"(a34)",end=9999)tmp
      if(tmp.ne.tag)goto 10

      idkh=1

9999  continue
      return
      end

c-----------------------------------------------------------------------
c--- check density fitting
c-----------------------------------------------------------------------
      subroutine rddf(idf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*16 tmp,tag

      imop=45
      tag=' Coulomb fitting'

      idf=0
      rewind(imop)
10    read(imop,"(a16)",end=9999)tmp
      if(tmp.ne.tag)goto 10

      idf=1

9999  continue
      return
      end

c-----------------------------------------------------------------------
c--- search dummy atoms
c-----------------------------------------------------------------------
      subroutine rddummy(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*14 tmp,tag

      imop=45
      tag=' Dummy atoms: '

      ist=0
      rewind(imop)
10    read(imop,"(a14)",end=9999)tmp
      if(tmp.ne.tag)goto 10

      ist=1

9999  continue
      return
      end

c-----------------------------------------------------------------------
c--- print out the upper triangular matrix
c-----------------------------------------------------------------------
      subroutine wrmat(mlst)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000,maxtri=2001000)   ! maxtri = maxbas * (maxbas + 1) / 2
      dimension trimat(maxtri),mlst(*)

c--- note: nbas=netot
      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)

      i47f=55

      write(i47f,"(' $OVERLAP')")
c--- DEBUG
c      if(mlst(4).ne.0) call romat(ovlp)
c--- reorder
      call romat(ovlp)
      call m2tri(netot,ovlp,trimat)
      call wrmat2(netottri,trimat)
      write(i47f,"(' $END')")

      write(i47f,"(' $DENSITY')")
c--- DEBUG
c      if(mlst(4).ne.0) call romat(dens)
c--- reorder
      call romat(dens)
      call m2tri(netot,dens,trimat) 
      call wrmat2(netottri,trimat)
      write(i47f,"(' $END')")

      if(mlst(3).eq.1)then
        write(i47f,"(' $FOCK')")
c--- DEBUG
c        if(mlst(4).ne.0) call romat(fock)
c--- reorder
        call romat(fock)
        call m2tri(netot,fock,trimat)
        call wrmat2(netottri,trimat)
        write(i47f,"(' $END')")
      end if

c--- MO is always square!!!
c--- orb(i,j); i for MO; j for AO basis; e.g. orb(1,j) = MO-1
      if(mlst(5).eq.1)then
        write(i47f,"(' $LCAOMO')")
c--- DEBUG
c        if(mlst(4).ne.0) call romat(orb)
c--- reorder
        call romat(orb)
        call wrsqmat(netot,orb)
        write(i47f,"(' $END')")
      end if

      if(mlst(6).eq.1)then
        write(i47f,"(' $KINETIC')")
c--- DEBUG
c        if(mlst(4).ne.0) call romat(ekin)
c--- reorder
        call romat(ekin)
        call m2tri(netot,ekin,trimat)
        call wrmat2(netottri,trimat)
        write(i47f,"(' $END')")
      end if

      if(mlst(7).eq.1)then
        write(i47f,"(' $NUCLEAR')")
c--- DEBUG
c        if(mlst(4).ne.0) call romat(epot)
c--- reorder
        call romat(epot)
        call m2tri(netot,epot,trimat)
        call wrmat2(netottri,trimat)
        write(i47f,"(' $END')")
      end if

      if(mlst(8).eq.1)then
        write(i47f,"(' $DIPOLE')")
c--- DEBUG
c        if(mlst(4).ne.0) call romat(dipx)
c--- reorder
c--- dip-x
        call romat(dipx)
        call m2tri(netot,dipx,trimat)
        call wrmat2(netottri,trimat)
c--- DEBUG
c        if(mlst(4).ne.0) call romat(dipy)
c--- reorder
c--- dip-y
        call romat(dipy)
        call m2tri(netot,dipy,trimat)
        call wrmat2(netottri,trimat)
c--- DEBUG
c        if(mlst(4).ne.0) call romat(dipz)
c--- reorder
c--- dip-z
        call romat(dipz)
        call m2tri(netot,dipz,trimat)
        call wrmat2(netottri,trimat)
        write(i47f,"(' $END')")
      end if

      return
      end

c-----------------------------------------------------------------------
c--- reorder A-matrix.
c-----------------------------------------------------------------------
      subroutine romat(a)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      common/mblock/nblk,nele(8),netot,netottri
      common/RELATION/lmap(maxbas)
      dimension a(maxbas,*),tmp(maxbas,maxbas)

      do i=1,netot
        j=lmap(i)
        tmp(i,1:netot)=a(j,1:netot)
      end do

      do i=1,netot
        a(i,1:netot)=tmp(i,1:netot)
      end do

      do i=1,netot
        j=lmap(i)
        tmp(1:netot,i)=a(1:netot,j)
      end do

      do i=1,netot
        a(i,1:netot)=tmp(i,1:netot)
      end do

      return
      end

c-----------------------------------------------------------------------
c--- print out the upper triangular matrix (2)
c-----------------------------------------------------------------------
      subroutine wrmat2(n,trimat)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension trimat(*)

      i47f=55

      write(i47f,"(2x,5e15.7e2)")(trimat(i),i=1,n)

      return
      end

c-----------------------------------------------------------------------
c--- write out the squere matrix
c-----------------------------------------------------------------------
      subroutine wrsqmat(n,a)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      dimension a(maxbas,maxbas),b(5)

      i47f=55

      line=1
      idx=0
100   do ip=1,5
        idx=idx+1
        if(idx.gt.n)then
          idx=idx-n
          line=line+1
          if(line.gt.n)goto 200
        end if
        b(ip)=a(line,idx)
      end do
      write(i47f,500)(b(i),i=1,5)
      goto 100

200   ip=ip-1
      if(ip.gt.0)write(i47f,500)(b(i),i=1,ip)

500   format(2x,5e15.7e2)
      return
      end

c-----------------------------------------------------------------------
c--- read Fock matrix
c-----------------------------------------------------------------------
      subroutine rdfock(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*200 tmp1
      character*26 tag1
      character*20 tag2

      imop=45
      tag1=' Closed-shell fock matrix '
      tag2=' made using density '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a200)",end=9999)tmp1
      itag=index(tmp1,tag1)*index(tmp1,tag2)
      if(itag.eq.0)goto 10

      ist=-1
      read(imop,"(a200)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,fock,imat)
      if(imat.eq.0)goto 9999

      ist=1

9999  return
      end

c-----------------------------------------------------------------------
c--- read natural orbital matrix, which is generated from density matrix.
c-----------------------------------------------------------------------
      subroutine rdmo2(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*200 tmp1
      character*18 tag1,tag4
      character*13 tag2
      character*9 tag3

      imop=45
      tag1=' Natural orbitals '
      tag2=' for density '
      tag3='# MATRIX '
      tag4=' ORBITALS NATURAL '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a200)",end=9999)tmp1
      itag=index(tmp1,tag1)*index(tmp1,tag2)
      if(itag.eq.0)goto 10
20    read(imop,"(a200)",end=9999)tmp1
      itag=index(tmp1,tag3)*index(tmp1,tag4)
      if(itag.eq.0)goto 20

      ist=-1
      call rdmat(nblk,nele,netot,orb,imat)
      if(imat.eq.0)goto 9999

      ist=1

9999  return
      end

c-----------------------------------------------------------------------
c--- read MO matrix
c-----------------------------------------------------------------------
      subroutine rdmo(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*200 tmp1
      character*10 tag1
      character*18 tag2

      imop=45
      tag1=' Orbitals '
      tag2=' read from record '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a200)",end=9999)tmp1
      if(tmp1(1:10).ne.tag1 .or. index(tmp1(11:200),tag2).eq.0)goto 10

      ist=-1
      read(imop,"(a200)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,orb,imat)
      if(imat.eq.0)goto 9999

      ist=1

9999  return
      end

c-----------------------------------------------------------------------
c--- read density matrix
c-----------------------------------------------------------------------
      subroutine rddens(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*200 tmp1
      character*9 tag1
      character*18 tag2

      imop=45
      tag1=' Density '
      tag2=' read from record '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a200)",end=9999)tmp1
      itag=index(tmp1,tag1)*index(tmp1,tag2)
      if(itag.eq.0)goto 10

      ist=-1
      read(imop,"(a200)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,dens,imat)
      if(imat.eq.0)goto 9999

      ist=1

9999  return
      end

c-----------------------------------------------------------------------
c--- read overlap matrix
c-----------------------------------------------------------------------
      subroutine rdovlp(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*26 tag,tmp1

      imop=45
      tag=' Overlap matrix loaded to '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a26)",end=9999)tmp1
      if(tmp1.ne.tag)goto 10

      ist=-1
      read(imop,"(a26)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,ovlp,imat)
      if(imat.eq.0)goto 9999

      ist=1
c--- check ovlp(i,i)
      call chkom(ist,netot,ovlp)
      if((ist.eq.2) .and. (nblk.eq.1))ifist=3    ! for C1 point group

9999  return
      end

c-----------------------------------------------------------------------
c--- read potential matrix
c-----------------------------------------------------------------------
      subroutine rdpote(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*22 tag1,tmp1

      imop=45
      tag1=' Operator V loaded to '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a22)",end=9999)tmp1
      if(tag1.ne.tmp1)goto 10

      ist=-1
      read(imop,"(a22)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,epot,imat)
      if(imat.eq.0)goto 9999

      ist=1

9999  return
      end

c-----------------------------------------------------------------------
c--- read kinetic matrix
c-----------------------------------------------------------------------
      subroutine rdkine(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*22 tag1,tmp1

      imop=45
      tag1=' Operator T loaded to '

      ist=0
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999
10    read(imop,"(a22)",end=9999)tmp1
      if(tag1.ne.tmp1)goto 10

      ist=-1
      read(imop,"(a22)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,ekin,imat)
      if(imat.eq.0)goto 9999

      ist=1

9999  return
      end

c-----------------------------------------------------------------------
c--- read dipole matrix (x,y,z)
c-----------------------------------------------------------------------
      subroutine rddip(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000,tol=1.d-15,con=0.529177249d0)  ! this con is used in NBO3

      common/mblock/nblk,nele(8),netot,netottri
      common/matrix/ovlp(maxbas,maxbas),dens(maxbas,maxbas),
     *fock(maxbas,maxbas),orb(maxbas,maxbas),ekin(maxbas,maxbas),
     *epot(maxbas,maxbas),dipx(maxbas,maxbas),dipy(maxbas,maxbas),
     *dipz(maxbas,maxbas)
      character*23 tag1,tmp1

      imop=45
      tag1=' Operator DM loaded to '

      ist=0
c--- symmetry can not be used
      if(nblk.ne.1)goto 9999

      ist=-1
      rewind(imop)
      call locmatrop(imat)
      if(imat.eq.0)goto 9999

c--- dip-x
10    read(imop,"(a23)",end=9999)tmp1
      if(tag1.ne.tmp1)goto 10

      ist=-2
      read(imop,"(a23)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,dipx,imat)
      if(imat.eq.0)goto 9999

      ist=-3
c--- dip-y
20    read(imop,"(a23)",end=9999)tmp1
      if(tag1.ne.tmp1)goto 20

      ist=-4
      read(imop,"(a23)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,dipy,imat)
      if(imat.eq.0)goto 9999

      ist=-5
c--- dip-z
30    read(imop,"(a23)",end=9999)tmp1
      if(tag1.ne.tmp1)goto 30

      ist=-6
      read(imop,"(a23)")tmp1
      if(tmp1(1:9).ne.'# MATRIX ')goto 9999
      call rdmat(nblk,nele,netot,dipz,imat)
      if(imat.eq.0)goto 9999

      ist=1
c--- scale
      do i=1,netot
        do j=1,netot
          dipx(j,i)=-dipx(j,i)*con
          if(abs(dipx(j,i)).lt.tol)dipx(j,i)=0.d0
          dipy(j,i)=-dipy(j,i)*con
          if(abs(dipy(j,i)).lt.tol)dipy(j,i)=0.d0
          dipz(j,i)=-dipz(j,i)*con
          if(abs(dipz(j,i)).lt.tol)dipz(j,i)=0.d0
        end do
      end do

9999  return
      end

c-----------------------------------------------------------------------
c--- check the overlap matrix
c-----------------------------------------------------------------------
      subroutine chkom(ist,netot,ovlp)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000,tol=1.d-9)
      dimension ovlp(maxbas,*)

      do i=1,netot
        if(abs(ovlp(i,i)-1.d0).gt.tol)then
          ist=2
          goto 100
        end if
      end do

100   return
      end

c-----------------------------------------------------------------------
c--- dump a matrix DATAMX into a upper triangular matrix DATATRI
c-----------------------------------------------------------------------
      subroutine m2tri(n,datamx,datatri)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      dimension datamx(maxbas,*),datatri(*)

      idx=0
      do i=1,n
        do j=1,i
          idx=idx+1
          datatri(idx)=datamx(j,i)
        end do
      end do
c--- DEBUG
c      write(8,"(10f10.5)")datatri(1:idx)
c---

      return
      end

c-----------------------------------------------------------------------
c--- read matrix
c--- imat=1 (overlap), 2 (density), 3 (fock)
c-----------------------------------------------------------------------
      subroutine rdmat(n,ne,nbas,datamx,imat)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      dimension ne(*),datamx(maxbas,*)
      
      imop=45
      do i=1,nbas
        datamx(1:nbas,i)=0.d0
      end do

      imat=1
      idx=0
      do i=1,n
        do j=1,ne(i)
          read(imop,"(5(f15.8,1x))",err=99)
     *    (datamx(k,idx+j),k=1+idx,ne(i)+idx)
        end do
        idx=idx+ne(i)
      end do
c--- DEBUG
c      do i=1,nbas
c        write(9,"(1000f10.5)")datamx(1:nbas,i)
c      end do
c---
      return

99    continue
      write(*,"(
     *' Error: format of the matrix is wrong!',/,
     *' Re-perform the MOLPRO calculation may solve this problem.',/)")
      imat=0
      return
      end

c-----------------------------------------------------------------------
c--- # b.s. in each block
c-----------------------------------------------------------------------
      subroutine rdblk(nbas,ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/mblock/nblk,nele(8),netot,netottri
      character*24 tag,tmp1

      imop=45
      tag=' NUMBER OF SYMMETRY AOS:'
      ist=1

      rewind(imop)
10    read(imop,"(a24)")tmp1
      if(tmp1.ne.tag)goto 10

      read(imop,"(41x,8(i4,5x))")(nele(i),i=1,nblk)
      netot=0
      do i=1,nblk
        netot=netot+nele(i)
      end do
      netottri=netot*(1+netot)
      netottri=netottri/2
c--- note: nbas=netot
      if(nbas.ne.netot)ist=0

      return
      end

c-----------------------------------------------------------------------
c--- read basis function types
c--- ist = -1 : can not find tag
c---        0 : normal
c---        1 : h,i,j... functions are found
c---        2 : symmetry equivalent atoms
c---  ig =  0 : no g-functions
c---        1 : g-functions are found
c-----------------------------------------------------------------------
      subroutine rdbstyp(ist,nbas,ig)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100,maxbas=2000)
      common/BASIS/icent(maxbas),label(maxbas)
      common/COORD/xyz(3,maxnatm),iz(2,maxnatm),natm
      character*22 tmp,tag
      character*5 clb(10)
      dimension ic(10)

      imop=45
      idx=0
      ig=0
      tag=' Basis function codes:'

      ist=-1
      rewind(imop)
10    read(imop,"(a22)",end=9999)tmp
      if(tmp.ne.tag)goto 10

      ist=0
      call lines(nbas,10,nline,last)

      do i=1,nline
        read(imop,"(1x,10(i3,1x,a5))")
     *ic(1),clb(1),ic(2),clb(2),ic(3),clb(3),ic(4),clb(4),ic(5),clb(5),
     *ic(6),clb(6),ic(7),clb(7),ic(8),clb(8),ic(9),clb(9),ic(10),clb(10)
        nj=10
        if(i.eq.nline)nj=last
        do j=1,nj
          icent(idx+j)=ic(j)
          label(idx+j)=maplab(clb(j))
c--- g-functions are found!
          if(label(idx+j).ge.400)ig=1
          if(label(idx+j).eq.0)then
c--- h or higher functions are found!
            write(*,"(' LABEL=',a5)")clb(j)
            ist=1
            goto 9999
          end if
        end do
        idx=idx+10
      end do

c--- check: if there are symmetry equivalent atoms, some of the atoms will be
c--- abcent in the atom list.
      idx=0
      do i=1,natm
        do j=1,nbas
          if(icent(j).eq.i)goto 100
        end do
        ist=2
100     if(ist.ne.0)goto 9999
      end do

9999  continue

      return
      end

c-----------------------------------------------------------------------
c--- the number of symmetry blocks
c-----------------------------------------------------------------------
      subroutine rdsymm(ist)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*17 tmp
      character*14 tag
      character*3 pglib(8),pg
      dimension nblib(8)
      data pglib/'C1 ','CS ','C2 ','CI ','D2 ','C2V','C2H','D2H'/
      data nblib/ 1,    2,    2,    2,    4,    4,    4,    8/
      common/mblock/nblk,nele(8),netot,netottri

      imop=45
      ist=0
      tag=' Point group  '

      rewind(imop)
10    read(imop,"(a17)",end=9999)tmp
      if(tmp(1:14).ne.tag)goto 10

100   ist=1
      pg=tmp(15:17)
      call charl2u(pg,3)
      do i=1,8
        if(pg.eq.pglib(i))then
          nblk=nblib(i)
          goto 9999
        end if
      end do

9999  continue

      return
      end

c-----------------------------------------------------------------------
c--- count the number of lines
c-----------------------------------------------------------------------
      subroutine lines(n1,n2,nline,last)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      nline=n1/n2
      last=n2
      if(mod(n1,n2).ne.0)then
        last=mod(n1,n2)
        nline=nline+1
      end if

      return
      end

c-----------------------------------------------------------------------
c--- map relationship AO type --> label
c--- 0: unknown or higher functions than g
c-----------------------------------------------------------------------
      function maplab(str)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(nty=60)
      character*5 str,typlib(nty)
      data typlib/
     *'s    ','x    ','y    ','z    ',                                    !  4 * cart s,p
     *'1s   ','2px  ','2py  ','2pz  ',                                    !  4 * pure s,p
     *'xx   ','xy   ','xz   ','yy   ','yz   ','zz   ',                    !  6 * cart d
     *'3d2- ','3d1+ ','3d1- ','3d2+ ','3d0  ',                            !  5 * pure d
     *'xxx  ','xxy  ','xxz  ','xyy  ','xyz  ','xzz  ','yyy  ','yyz  ',    ! 10 * cart f
     *'yzz  ','zzz  ',                                                    !
     *'4f0  ','4f1+ ','4f1- ','4f2+ ','4f2- ','4f3+ ','4f3- ',            !  7 * pure f
     *'xxxx ','xxxy ','xxxz ','xxyy ','xxyz ','xxzz ','xyyy ','xyyz ',    ! 15 * cart g
     *'xyzz ','xzzz ','yyyy ','yyyz ','yyzz ','yzzz ','zzzz ',            !
     *'5g0  ','5g1+ ','5g1- ','5g2+ ','5g2- ','5g3+ ','5g3- ','5g4+ ',    !  9 * pure g
     *'5g4- '/
      dimension idx(nty)
      data idx/
     *001,101,102,103,
     *001,101,102,103,
     *201,202,203,204,205,206,
     *251,252,253,254,255,
     *301,302,303,304,305,306,307,308,309,310,
     *351,352,353,354,355,356,357,
     *401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,
     *451,452,453,454,455,456,457,458,459/

      it=0

      do i=1,nty
        if(str.eq.typlib(i))then
          it=idx(i)
          goto 9999
        end if
      end do

9999  continue
      maplab=it

      return
      end

c-----------------------------------------------------------------------
c--- read #CONTRACTION
c--- nbas < 0: #groups of basis
c--- nbas > 0: #basis
c-----------------------------------------------------------------------
      subroutine rdnbas(nbas,ist,idf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      character*24 tmp,tag

      imop=45
      tag=' NUMBER OF SYMMETRY AOS:'
      ist=0

      nbas=0
      rewind(imop)
10    read(imop,"(a24)",end=99)tmp
      if(tmp.eq.tag) nbas=nbas-1
      goto 10
99    if( (nbas.ne.-1 .and. idf.eq.0) .or.
     *    (nbas.lt.-2 .and. idf.eq.1) .or.
     *    (nbas.eq.0) ) goto 9999

      rewind(imop)
100   read(imop,"(a24)")tmp
      if(tmp.eq.tag) goto 200
      goto 100

200   read(imop,"(30x,i7)")nbas
      if(nbas.gt.maxbas)then
        write(*,"(/,' NBAS=',i8,',  MAXBAS=',i8)")nbas,maxbas
        ist=1
      end if

9999  continue
      return
      end

c-----------------------------------------------------------------------
c--- plot orbitals (1) or not (0)
c-----------------------------------------------------------------------
      subroutine plotnbo(ipl)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*1 yn

100   write(*,"(/,
     *' Do you want to plot orbitals generated by NBO?',/,
     *' [Y or N (default)]',/,
     *' > ',$)")
      read(*,"(a1)")yn
      lenth=len_trim(yn)
      if(lenth.eq.0.or.yn.eq.'N'.or.yn.eq.'n')then
        ipl=0
        write(*,"(//,' The PLOT keyword is OFF.')")
      else if(yn.eq.'Y'.or.yn.eq.'y')then
        ipl=1
        write(*,"(//,' The PLOT keyword is ON.')")
      else
        write(*,1000)yn
        goto 100
      end if

      return
1000  format(' Unknown "',a1,'". Try again.')
      end

c-----------------------------------------------------------------------
c--- head
c-----------------------------------------------------------------------
      subroutine head(ver,dt)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*10 dt
      character*5 ver

      write(*,"(//,1x,45('*'),/
     *'           * * *     MOLBO     * * *',/,
     *'          Version ',a5,',   ',a10,//,
     *'    An interface of MOLPRO to NBO and MBO.',/,
     *1x,45('*'),/)")ver,dt

      return
      end

c-----------------------------------------------------------------------
c--- datalist
c-----------------------------------------------------------------------
      subroutine datalist
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(nop=4)
      character*7 oplist(nop)
      data oplist/'CORE   ','CHOOSE ','DEL    ','NRTSTR '/

      write(*,"(//,
     *' You can also define the following datalists,')")
      do i=1,nop
        if(i.ne.3)then
          write(*,100)i,oplist(i)
        else
          write(*,200)i,oplist(i)
        end if
      end do

      write(*,300)

      return

100   format(2x,i1,') $',a7)
200   format(2x,i1,') $',a7,' (RHF/UHF only; see $FOCK)')
300   format(//,
     *" *** NOTE ***",//,
     *" To reproduce the NBO results by other programs, the exactly",/,
     *" same method and coordinate orientation should be used.")
      end

c-----------------------------------------------------------------------
c--- write $BASIS
c-----------------------------------------------------------------------
      subroutine wrbas(nbas)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxbas=2000)
      common/BASIS/icent(maxbas),label(maxbas)
      character*10 tag(3)
      data tag/'  CENTER =','   LABEL =','          '/

      i47f=55

      write(i47f,"(' $BASIS')")

      nstep=17
      call lines(nbas,nstep,nline,last)
c--- CENTER
      call iwrite(0,tag(1),tag(3),i47f,nline,last,nstep,icent)
c--- LABEL
      call iwrite(0,tag(2),tag(3),i47f,nline,last,nstep,label)

      write(i47f,"(' $END')")

1110  format(a10,17i4)

      return
      end

c-----------------------------------------------------------------------
c--- write NBO keywords
c-----------------------------------------------------------------------
      subroutine wrnbokey(nbas,iplot)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100)
      common/COORD/xyz(3,maxnatm),iz(2,maxnatm),natm
      parameter(maxexp=5000,maxgroup=800)
      common/CONTRACT/NSHELL,NEXP,NCOMP(maxgroup),NPRIM(maxgroup),   ! for NBO
     *NPTR(maxgroup),EX(maxexp),CS(maxexp),CP(maxexp),CD(maxexp),
     *CF(maxexp),CG(maxexp)
      character*4 cnatm
      character*5 cnbas

      i47f=55
      i31f=56

      write(cnatm,"(i4)")natm
      cnatm=ADJUSTL(cnatm)
      write(cnbas,"(i5)")nbas
      cnbas=ADJUSTL(cnbas)

c--- .47
      write(i47f,"(' $GENNBO  NATOMS=',a4,' NBAS=',a5,
     *' UPPER  BODM  $END')")cnatm,cnbas
      if(iplot.eq.0)then
        write(i47f,"(' $NBO BNDIDX NLMO $END')")
      else
        write(i47f,"(' $NBO BNDIDX NLMO PLOT $END')")
      end if

c--- .31
      if(iplot.lt.0)then
        write(i31f,"(
     *' Use this file to overwrite the *.31 file generated by NBO 3.',/,
     *' Basis set information needed for plotting orbitals',/,
     *1x,75('-'))")
        write(i31f,"(1x,3i6,/,1x,75('-'))")natm,NSHELL,NEXP
      end if

      return
      end

c-----------------------------------------------------------------------
c--- write coordinates
c-----------------------------------------------------------------------
      subroutine wrcoord(ig,iplot)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100)
      common/COORD/xyz(3,maxnatm),iz(2,maxnatm),natm

      i47f=55
      i31f=56

c--- .47
      write(i47f,"(' $COORD')")
      if(ig.ne.0 .and. iplot.gt.0)then
        write(i47f,"(
     *' For NBO 4 and higher; generated by MOLBO.')")
      else
        write(i47f,"(
     *' For NBO 3 and higher; generated by MOLBO.')")
      end if
c     About NBO6:
c     It may lead to numerical errors of about 1.0d-6 in the overlap matrix,
c     which cannot pass the examination of NBO6. More digits should be printed.
      do i=1,natm
cooo        write(i47f,"(1x,2i5,3f15.6)")iz(:,i),xyz(:,i)
        write(i47f,"(1x,2i5,3f18.9)")iz(:,i),xyz(:,i)
      end do
      write(i47f,"(' $END')")

c--- .31
      if(iplot.lt.0)then
        do i=1,natm
          write(i31f,"(i5,3f14.9)")iz(1,i),xyz(:,i)
        end do
        write(i31f,"(1x,75('-'))")
      end if

      return
      end

c-----------------------------------------------------------------------
c--- read coordinates
c-----------------------------------------------------------------------
      subroutine rdcoord(ird,idf)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnatm=100)
      common/COORD/xyz(3,maxnatm),iz(2,maxnatm),natm
      character*19 tmp
      character*2 elm

      imop=45
      c=0.529177249d0  ! this con is used in NBO3

      ird=0
      rewind(imop)
10    read(imop,"(19a)",end=99)tmp
      if(tmp.eq.' ATOMIC COORDINATES') ird=ird+1
      goto 10
99    if( (ird.ne.1 .and. idf.eq.0) .or.
     *    (ird.gt.2 .and. idf.eq.1) )goto 9999

      rewind(imop)
100   read(imop,"(19a)",end=9999)tmp
      if(tmp.eq.' ATOMIC COORDINATES') goto 200
      goto 100

200   read(imop,"(//)",end=9999)

240   read(imop,"(i4,2x,a2,2x,f8.2,3f15.9)",end=9999)
     *  iatm,elm,zatm,xyz(:,iatm)
      if(iatm.gt.maxnatm)goto 9990
      if(iatm.eq.0)goto 9999

      call charl2u(elm,2)
      iz(1,iatm)=izatm(elm)
      iz(2,iatm)=int(zatm)
c--- Bohr --> Ang. The unit of this part in molpro is always Bohr.
      xyz(:,iatm)=xyz(:,iatm)*c
      natm=iatm
      goto 240

9990  continue
      write(*,"(' IATM=',i4,',  MAXNATM=',i4)")iatm,maxnatm
      ird=-1
9999  continue
c--- DEBUG
c      do i=1,natm
c        write(*,"(2i5,3f15.6)")iz(:,i),xyz(:,i)
c      end do
c
      return
      end

c-----------------------------------------------------------------------
c--- return IZ of element EL
c-----------------------------------------------------------------------
      function izatm(EL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*2 EL
      parameter (max_za=103)
      character*2 ATOMLIB(max_za)
      data (atomlib(i),i=1,50) /
     1'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     2'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     3'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     4'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     5'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN'/
      data (atomlib(i),i=51,100) /
     1'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
     2'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     3'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
     4'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     5'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM'/
      data (atomlib(i),i=101,103) /
     1'MD','NO','LR'/

c--- 1. EL(2:2) may be a number. Replace it by a space.
      ic=ichar(EL(2:2))
      if(ic.ge.48.and.ic.le.57)EL(2:2)=' '

c--- 2. get IZ
      do i=1,max_za
        if(EL.eq.atomlib(i))goto 100
      end do
      write(*,99)EL
      stop
99    format(" Wrong! Unknown element '",a2,"'")

100   izatm=i

      return
      end

c-----------------------------------------------------------------------
c--- define MOLPRO output, *.31, *.47, and mbo files
c-----------------------------------------------------------------------
      subroutine molout(f31,f47,fmb)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*27 f31,f47,fmop(2)
      character*7 exten(4)
      character*32 fmb
      data exten /'.out   ','.OUT   ','.output','.OUTPUT'/

      imop=45

      write(*,"(/)")
100   write(*,"(' MOLPRO output file name within 20 characters:')")
      write(*,"(' > ',$)")
      read(*,"(a20)")fmop(1)(:)
      lstr=nonspace(fmop(1)(:))
      lend=LEN_TRIM(fmop(1)(:))
      if(lend.eq.0)then                 ! use default file name
        lstr=1
        lend=6
        fmop(1)(1:6)='MOLPRO'
      end if
      open(imop,file=fmop(1)(lstr:lend),status='old',err=110)
      iinp=1
      goto 300
110   if(fmop(1)(lend:lend).eq.'.')lend=lend-1
      iinp=2
      do i=1,4
        fmop(2)(:)=fmop(1)(lstr:lend)//trim(exten(i))
        open(imop,file=fmop(2)(:),status='old',err=120)
        goto 300
120     continue
      end do
      write(*,"(//,' These MOLPRO output files do not exist!')")
      write(*,"(2x,a)")fmop(1)
      do i=1,4
        write(*,"(2x,a)")fmop(1)(lstr:lend)//trim(exten(i))
      end do
      write(*,"(/,' Please try again.',/)")
      goto 100

300   write(*,"(/,' The MOLPRO output ',a,' has been found.',/,
     *1x,45('-'),//)")
     *trim(fmop(iinp))

c--- define the *.31 file name
      f31=fmop(iinp)(lstr:lend)//'.31'
c--- define the *.47 file name
      f47=fmop(iinp)(lstr:lend)//'.47'
c--- define the *.mbo file name
      fmb=fmop(iinp)(lstr:lend)//'-mbo.out'

      return
      end

c-----------------------------------------------------------------------
c--- locate ' PROGRAM * MATROP'
c-----------------------------------------------------------------------
      subroutine locmatrop(ifd)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*17 tag,tmp

      imop=45
      ifd=0
      tag=' PROGRAM * MATROP'
1     read(imop,"(a17)",end=9999)tmp
      if(tmp.ne.tag)goto 1
      ifd=1

9999  return
      end

c-----------------------------------------------------------------------
c--- tmp --> TMP
c-----------------------------------------------------------------------
      subroutine charl2u(tmp,nc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*1 tmp(*),L2U

      do i=1,nc
        tmp(i)=L2U(tmp(i))
      end do

      return
      end

c-----------------------------------------------------------------------
c--- l-->L
c-----------------------------------------------------------------------
      function L2U(letter)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*1 letter,L2U
      if((ichar(letter).ge.97).and.(ichar(letter).le.122))then
        L2U=char(ichar(letter)-32)
      else
        L2U=letter
      endif
      return
      end

c-----------------------------------------------------------------------
c--- length of a string without the first and last spaces.
c-----------------------------------------------------------------------
      function nonspace(string)
      implicit double precision (a-h,o-z)
      character*(*) string
      character*1 space

      space=' '
      length=LEN_TRIM(string)
      if(length.eq.0) then
       i=1
      else
       do i=1,length
         if(string(i:i).ne.space) goto 20
       end do
      endif

20    nonspace=i

      return
      end

c-----------------------------------------------------------------------
c--- read an <ENTER> and then stop
c-----------------------------------------------------------------------
      subroutine estop
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write(*,"(//,' Press <ENTER> to exit',/)")
      read(*,*)

      stop

      return
      end

c-----------------------------------------------------------------------
