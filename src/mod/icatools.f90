module icatools
  ! Author: Davide Maino
  ! CHANGES:
  ! original version 21/MAR/2001 by D.Maino@OAT
  ! add remove_dipole_custom: 19 Jan 2004 DMM@UniMi
  use healpix_types
  use pix_tools, ONLY: nside2npix, pix2vec_ring, pix2vec_nest
contains
  !========================================================================
  subroutine remove_dipole_custom(nside, map, nmaps, ordering, degree, &
       & multipoles, cut_pixel, fmissval)
    !======================================================================
    ! removes monopole (and dipole) from a map

    !
    ! Nside:     I4,       IN   : Healpix resolution parameter
    ! map:       R4, array,INOUT: Heapix map (see Notes below)
    ! ordering:  I4,       IN:   Healpix scheme 1:RING, 2: NESTED
    ! degree:    I4,       IN:   multipole to remove, 1: monopole, 2: monopole and dipole
    ! multipoles:R8, array,OUT:  value of monopole and dipole
    ! fmissval:  R4, Option, IN: value used to flag bad pixel on input, default=-1.6375e30
    !
    ! note : if degree= 1, or 2, the map is modified on output
    !     * the monopole (and dipole) is/are removed
    !     * pixels within the symmetric cut parameterized 
    !       by cos_theta_cut are set to fmissval (or its default value)
    !  if degree = 0, nothing is done
    !  all other values of degree are invalid
    !
    ! v1.0, EH, Caltech, Jan-2002, based on homonyme IDL routine
    ! v2.0, DM, UniMi, Jan-2004, modified for custom galactic cuts
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    implicit none
    integer(kind=i4b),                  intent(in)    :: nside
    real   (kind=sp), dimension(0:),    intent(inout) :: map
    integer(kind=i4b),                  intent(in)    :: nmaps
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP), dimension(0:),    intent(out)   :: multipoles
    real   (kind=SP),                   intent(in), optional :: fmissval


    integer(kind=i4b)                 :: ipix, npix, n_mult
    integer(kind=i4b)                 :: i
    logical(lgt)                      :: dodipole
    real(kind=dp)                     :: flag, temp, wmin
    real(kind=dp), dimension(1:3)     :: vec
    real(kind=dp), dimension(0:3)     :: b, wdiag
    real(kind=dp), dimension(0:3,0:3) :: mat, umat, vmat
    real(kind=dp), dimension(0:3)     :: dmultipoles
    real(kind=SP)                     :: fmiss_effct
    character(len=*), parameter :: code = "REMOVE_DIPOLE"
    real(kind=sp),    parameter :: fbad_value = -1.6375e30_sp
!    real(kind=dp)                     :: theta1, theta2
!    integer(kind=i4b)                 :: ncpix, ncp
    integer(kind=i4b), dimension(1:)  :: cut_pixel
    !============================================================

    npix = nside2npix(nside)
    multipoles = 0.0_dp
    n_mult = (degree)**2
    fmiss_effct = fbad_value
    if (present(fmissval)) fmiss_effct = fmissval

    if (degree == 0) then
       print*," No monopole nor dipole removal"
       return
    elseif (degree == 1) then
       dodipole = .false.
    else if (degree == 2) then
       dodipole = .true.
    else
       print*,code//"> degree can only be "
       print*,"      1: monopole (l=0) removal or "
       print*,"      2: monopole and dipole (l=0,1) removal"
       print*,code//"> ABORT ! "
       stop
    endif

    !----------------------------------------------
    ! flag out pixels within specified cut
    !----------------------------------------------
    map(cut_pixel(:)) = fmiss_effct
    
    !----------------------------------------------
    ! generate least square linear system
    !----------------------------------------------
    mat = 0.0_dp
    b   = 0.0_dp
    do ipix = 0, npix-1

       ! flag = 1 for good values
       !      = 0 for bad values
       flag = 1.0_dp
       !     if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) flag = 0.0_dp
       if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) goto 20

       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
       endif

       ! construct vector T*(1,x,y,z)
       temp = map(ipix) * flag
       b(0) = b(0) + temp
       if (dodipole) then
          ! computes dipole basis functions
          b(1:3) = b(1:3) + temp * vec(1:3)
       endif

       ! construct matrix (1,x,y,z)#(1,x,y,z)
       mat(0,0) = mat(0,0) + flag
       if (dodipole) then
          do i = 1, 3
             mat(i,0)   = mat(i,0)   + vec(i) * flag
             mat(i,1:3) = mat(i,1:3) + vec(i) * vec(1:3) * flag
          enddo
       endif

20     continue
    enddo

    ! first row = first column (and vice versa)
    mat(0,1:3) = mat(1:3,0)


    !----------------------------------------------
    ! solve system    mat . (mono, dip_x, dip_y, dip_z) = b
    !----------------------------------------------

    if (dodipole) then
       ! SVD decomposition
       umat = mat
       call dsvdcmp(umat, 4, 4, 4, 4, wdiag, vmat)
       ! thresholding
       wmin = maxval(wdiag)* 1.e-6_dp
       where (wdiag < wmin)
          wdiag = 0.0_dp
       end where
       ! back substitution
       call dsvbksb(umat, wdiag, vmat, 4, 4, 4, 4, b, dmultipoles)

    else
       dmultipoles(0) = b(0) / mat(0,0) ! average temperature
    endif

    !----------------------------------------------
    ! remove monopole and dipole
    !----------------------------------------------
    do ipix = 0, npix-1

       if ( abs(map(ipix) - fmiss_effct) < 1.e-5*fmiss_effct ) goto 10

       map(ipix) = map(ipix) - dmultipoles(0)
       if (dodipole) then
          ! computes dipole basis functions
          ! pixel -> vector
          if (ordering == 1) call pix2vec_ring( nside, ipix, vec)
          if (ordering == 2) call pix2vec_nest( nside, ipix, vec)
          map(ipix) = map(ipix) - SUM( dmultipoles(1:3) * vec(1:3))
       endif

10     continue

    enddo

    multipoles(0:n_mult-1) = dmultipoles(0:n_mult-1)

    !

    return

  end subroutine remove_dipole_custom

  !=====================================================
  subroutine fpica(code, nsig, ncount, Xsig, S, whiteningMatrix, &
       & dewhiteningMatrix, g, a1, a2, epsilon, maxNumIterations, A, W, conv)
    !======================================================
    IMPLICIT NONE
    
    INTEGER(I4B) :: i,j,k,nsig,l,ncount, conv
    INTEGER(I4B) :: maxNumIterations,round,numOfIC
    INTEGER(I4B) :: numFailures,failureLimit,useNlinearity

    REAL(DP) :: a1,a2,epsilon,normdiff,normsum
    REAL(DP), DIMENSION(:),ALLOCATABLE :: wvet,wold,wdiff,wsum
    REAL(DP), DIMENSION(:),ALLOCATABLE :: usig,hypTan,gauss,dGauss
    REAL(DP), DIMENSION(nsig,ncount) :: Xsig
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: XsigT
    REAL(DP), DIMENSION(nsig,nsig) ::s,whiteningMatrix,dewhiteningMatrix,a,w
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: B,Ident,Bt
    CHARACTER(LEN=*) :: g, code

    AllOCATE(usig(ncount),B(nsig,nsig),Ident(nsig,nsig))
    ALLOCATE(wvet(nsig),wold(nsig),wdiff(nsig),wsum(nsig))
    ALLOCATE(Bt(nsig,nsig),xsigt(ncount,nsig))

    IF(g=='p')THEN
       useNlinearity=1
       PRINT*,"      "//code//"> Used nonlinearity [pow3]."   
    ELSEIF(g=='t') THEN
       useNlinearity=2
       PRINT*,"      "//code//"> Used nonlinearity [tanh]."
    ELSEIF(g=='l') THEN
       useNlinearity=4
       PRINT*,"      "//code//"> Used nonlinearity [pow4]."
    ELSEIF(g=='g') THEN
       useNlinearity=3
       PRINT*,"      "//code//"> Used nonlinearity [gauss]."
    ELSE
       PRINT*,"      "//code//"> ILLEGAL VALUE [ ',g,' ] for parameter: g."
    ENDIF
    PRINT*,"      "//code//"> Stopping criterion:",epsilon
    Ident=0.d0
    DO i=1,nsig
       DO j=1,nsig
          IF(i==j) Ident(i,j)=1.d0
       ENDDO
    ENDDO
    failureLimit=10
    numOfIC=nsig

    B=0.d0
    round=1
    numFailures=0
    COMP : DO i=1,nsig
       IF(round>numOfIC) EXIT
       call random_number(wvet)
!	print*,wvet
       wvet=wvet-0.5d0
       DO l=1,nsig
          DO k=1,nsig
             Bt(k,l) = B(l,k)
          ENDDO
       ENDDO
       wvet=wvet-MATMUL(B,MATMUL(Bt,wvet))
       wvet=wvet/sqrt(SUM(wvet(:)**2))
       wold=0.d0
       INDP : DO j=1,maxNumIterations+1
          conv = 1
          IF(j==maxNumIterations+1) THEN
             PRINT*,"      "//code//"> Component number ",round," did not converge in "&
                  &,maxNumIterations," iterations."
             conv = 0
             round=round-1
             numFailures=numFailures+1
             IF(numFailures>failureLimit) THEN
                PRINT*,"      "//code//"> Too many failures to conv.",numFailures,". Giving up."
                RETURN
             ENDIF
             EXIT INDP
          ENDIF
          DO l=1,nsig
             DO k=1,nsig
                Bt(k,l) = B(l,k)
             ENDDO
          ENDDO
          wvet=wvet-MATMUL(B,MATMUL(Bt,wvet))
          wvet=wvet/sqrt(SUM(wvet(:)**2))
          wdiff=wvet-wold
          wsum=wvet+wold
          normdiff=sqrt(SUM(wdiff(:)**2))
          normsum=sqrt(SUM(wsum(:)**2))
!	  print*,"normdiff = ",normdiff,"  normsum = ",normsum
          IF(normdiff<epsilon.or.normsum<epsilon) THEN
             numFailures=0
             DO k=1,nsig 
                B(k,round)=wvet(k)
                A(k,round)=SUM(dewhiteningMatrix(k,:)*wvet(:))
                W(round,k)=SUM(wvet(:)*whiteningMatrix(:,k))
             ENDDO
             PRINT*,"      "//code//"> IC [",i,"] computed in(",j,"steps)"
             EXIT INDP
          ENDIF
          wold=wvet
          DO l=1,nsig
             DO k=1,ncount
                XsigT(k,l) = Xsig(l,k)
             ENDDO
          ENDDO
          usig=MATMUL(XsigT,wvet)
          IF (useNlinearity==1) THEN 
             wvet=(MATMUL(Xsig,usig**3)-3.d0*SUM(usig(:)**2)*MATMUL(Ident+S,wvet))/ncount
          ELSEIF (useNlinearity==4) THEN
             wvet=(MATMUL(Xsig,usig**5)-5.d0*SUM(usig(:)**4)*MATMUL(Ident+S,wvet))/ncount
          ELSEIF(useNlinearity==2) THEN
             ALLOCATE(hypTan(ncount))
             hypTan=tanh(a1*usig)
             wvet=(MATMUL(Xsig,hypTan)-a1*SUM(hypTan(:))*MATMUL(Ident+S,wvet))/ncount
             DEALLOCATE(hypTan)
          ELSEIF(useNlinearity==3) THEN
             ALLOCATE(gauss(ncount),dGauss(ncount))
             DO k=1,ncount
                gauss(k)=usig(k)*EXP(-(a2/2.d0)*usig(k)**2)
                dGauss(k)=(1-a2*usig(k)**2)*EXP(-(a2/2.d0)*usig(k)**2)
             ENDDO
             wvet=(MATMUL(Xsig,gauss)-SUM(dGauss(:))*MATMUL(Ident+S,wvet))/ncount
             DEALLOCATE(gauss,dGauss)   
          ELSE
             PRINT*,"      "//code//"> Code for desidered nonlinearity not found!"
             RETURN
          ENDIF
          wvet=wvet/sqrt(SUM(wvet(:)*wvet(:)))
       ENDDO INDP
       round=round+1
     ENDDO COMP
     DEAllOCATE(usig,B,wvet,wold,wdiff,wsum,Ident)
     RETURN
  
  end subroutine fpica

  !===========================================================
  subroutine whitening(code,nsig,ncount,xsig,s,whiteningMatrix,dewhiteningMatrix)
    !===============================================================
    USE nrtype

    IMPLICIT NONE

    INTEGER(I4B) :: i,j,nsig,nrot,ncount

    REAL(DP),DIMENSION(:), ALLOCATABLE :: EINV
    REAL(DP),DIMENSION(nsig,ncount) :: xsig
    REAL(DP),DIMENSION(:,:), ALLOCATABLE :: xsigt,dummy
    REAL(DP),DIMENSION(nsig,nsig) :: s
    REAL(DP),DIMENSION(nsig,nsig) :: whiteningMatrix,dewhiteningMatrix
    REAL(DP),DIMENSION(:,:), ALLOCATABLE :: CORR,D,TEMP
    CHARACTER(LEN=*) :: code

    ALLOCATE(corr(nsig,nsig),d(nsig,nsig),temp(nsig,nsig))
    ALLOCATE(einv(nsig),xsigt(ncount,nsig),dummy(nsig,ncount))
!	print*,xsig(1,1:10)
    DO i=1,nsig
       DO j=1,ncount
          xsigt(j,i) = xsig(i,j)
       ENDDO
    ENDDO
    CORR=MATMUL(Xsig,Xsigt)
    CORR=CORR/ncount
!    print*,corr
    TEMP=CORR-S
!    print*,S
    call JACOBI(TEMP,EINV,CORR,nrot,nsig)
    do i=1,nsig
      if ((einv(i) < 0.) .and. (abs(einv(i)) < 1.e-7)) einv(i)=abs(einv(i))
    enddo
    PRINT*,"      "//code//"> These are the Principal Components:"
    PRINT*,EINV
    D=0.d0
    DO i=1,nsig
       DO j=1,nsig
          IF(EINV(i)<=0.d0) THEN
             PRINT*,"      "//code//">Negative eigenvalue of covariance matrix"
             STOP
          ENDIF
          IF(i==j) D(i,i)=sqrt(1.d0/EINV(i))
       ENDDO
    ENDDO
    whiteningMatrix=MATMUL(CORR,MATMUL(D,transpose(CORR)))

    dummy=MATMUL(whiteningMatrix,Xsig)
    Xsig = dummy

    D=0.d0
    DO i=1,nsig
       DO j=1,nsig
          IF(i==j) D(i,i)=sqrt(EINV(i))
       ENDDO
    ENDDO
    dewhiteningMatrix=MATMUL(MATMUL(CORR,D),transpose(CORR))

    S=MATMUL(MATMUL(whiteningMatrix,S),whiteningMatrix)

    DEALLOCATE(CORR,EINV,D,TEMP)
    RETURN
  end subroutine whitening

  !=========================================
  subroutine jacobi(a,d,v,nrot,nsig)
    !==========================================
    USE nrtype 
    USE nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,upper_triangle
    
    IMPLICIT NONE
    INTEGER(I4B)                           :: nrot
    REAL(DP), DIMENSION(nsig)    :: d
    REAL(DP), DIMENSION(nsig,nsig)  :: a
    REAL(DP), DIMENSION(nsig,nsig)  :: v
    INTEGER(I4B)                           :: i,ip,iq,n,nsig
    REAL(DP) :: c,g,h,s,sm,t,tau,theta,tresh
    REAL(DP), DIMENSION(size(d)) :: b,z
    n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')
    call unit_matrix(v(:,:))
    b(:)=get_diag(a(:,:))
    d(:)=b(:)
    z(:)=0.d0
    nrot=0
    do i=1,50
       sm=sum(abs(a),mask=upper_triangle(n,n))
       if (sm == 0.d0) RETURN
       tresh=merge(0.2d0*sm/n**2,0.d0, i < 4 )
       do ip=1,n-1
          do iq=ip+1,n
             g=100.d0*abs(a(ip,iq))
             if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                  .and. (abs(d(iq))+g == abs(d(iq)))) then
                a(ip,iq)=0.d0
             else if (abs(a(ip,iq)) > tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g == abs(h)) then
                   t=a(ip,iq)/h
                else
                   theta=0.5d0*h/a(ip,iq)
                   
                   t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                   if (theta < 0.d0) t=-t
                end if
                c=1.d0/sqrt(1+t**2)
                s=t*c
                tau=s/(1.d0+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.d0
                call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
                call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
                call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
                call jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             end if
          end do
       end do
       b(:)=b(:)+z(:)
       d(:)=b(:)
       z(:)=0.d0
     end do
     call nrerror('too many iterations in jacobi')
   CONTAINS
!BL
     SUBROUTINE jrotate(a1,a2)
       REAL(KIND(0.d0)), DIMENSION(:)        :: a1,a2
       REAL(KIND(0.d0)), DIMENSION(size(a1)) :: wk1
       wk1(:)=a1(:)
       a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
       a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
     END SUBROUTINE jrotate
  end subroutine jacobi

end module icatools





