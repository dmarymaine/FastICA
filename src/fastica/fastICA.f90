!--------------------------------------------------------
! This code is the main code for the fastICA analysis of 
! microwave sky maps. The hearth of the code is the 
! subroutines spherx, whitening and fpica.
! The first two make the input "sphered" in the sense
! that remove the mean from the input and make them
! white which means that make a linear transformation on
! the input to get their uncorrelated (white) with
! correlation matrix which is the identity matrix.
! This is also true in the case when you account for
! noise since it uses the "quasi"-whitening procedure.
! The FastICA algorithm is in subroutine fpica.
!
! The author of this code is Andrea Farusi (IEI-CNR Pisa)
! The present version has been modified by Davide Maino and
! Carlo Baccigalupi.
!
!--------------------------------------------------------
PROGRAM fastica
  !------------------------------------------------------
  ! EXTERNAL LIBRARIES:
  !    this code uses the HEALPix package 
  !      by Hivon and Gorski
  !
  !    this code uses also the CFITSIO library for
  !      handling FITS file
  !
  ! HISTORY:
  !    Jan 1999: first version: both QUAD-CUBE and HEALPIX
  !    Mar 2002: implemented for polarisation (Q & U separately)
  !              and frequency scaling reconstruction
  !    Dec 2002: implemented for polarisation (Q+U maps)
  !    Jan 2004: adapted for HEALPix_1.20
  !
  !------------------------------------------------------
  !    version 1.0.0
  !------------------------------------------------------
  USE healpix_types
  USE fitstools, ONLY : getsize_fits, input_map, read_par, read_dbintab, write_bintab
  USE head_fits, ONLY : add_card
  USE utilities, ONLY : die_alloc
  USE pix_tools, ONLY : npix2nside, vec2ang, remove_dipole
  USE extension, ONLY : getEnvironment, getArgument, nArguments
  USE paramfile_io, ONLY : paramfile_handle, parse_init, parse_int, &
       parse_string, parse_double, parse_lgt, concatnl
  USE icatools
  USE utils, ONLY : get_gen_cut

  IMPLICIT NONE
  INTEGER(I4B) :: nsmax, mlpol,ipix
  INTEGER(I4B) :: i, j, k, l, kk, m, ll
  INTEGER(I4B) :: nsig, npixtot, nmaps, ncount, ncut
  INTEGER(I4B) :: icatype, qutype
  INTEGER(I4B) :: nstokes, nsample
  INTEGER(I4B) :: iseed, iargc
  INTEGER(I4B) :: ordering, order_type
  INTEGER(I4B) :: nlheader, conv
  INTEGER(I4B) :: choice_rmsprior, choice_cut, cure_noise
  INTEGER(I4B) :: maxNumIterations
  INTEGER(I4B) :: status, type
  INTEGER(I4B), DIMENSION(8,2) :: values_time
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: cut_pixel, left_pixel

  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: mask  

  REAL(DP) :: a1,a2,epsilon, cos_theta_cut
  REAL(DP) :: b_min,b_max,l_min,l_max
  REAL(KIND=SP) :: fmissval
  REAL(KIND=SP) :: clock_time
  REAL(KIND=DP) :: nullval
  REAL(KIND=SP), PARAMETER :: fbad_value = -1.6375e30_sp

  REAL(DP), DIMENSION(:),   ALLOCATABLE :: cl, cl_E, cl_B, freq, med
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: rms_guess
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sig, s, xsig, n, errmap
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: wmat, dwmat, a, w
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: mapIO, err
  real(kind=DP), dimension(1:2) :: zbounds

  CHARACTER(LEN=filenamelen) :: g 
  CHARACTER(LEN=3) :: output
  CHARACTER(LEN=filenamelen), DIMENSION(:,:), ALLOCATABLE :: files, noises, rmsfiles
  CHARACTER(LEN=80), DIMENSION(1:120) :: header
  LOGICAL(LGT) :: oknoise, okrms, noise_in_data
  LOGICAL(LGT) :: polarisation

  CHARACTER(LEN=filenamelen) :: description
  CHARACTER(LEN=filenamelen) :: paramfile
  CHARACTER(LEN=filenamelen) :: maskfile 
  CHARACTER(LEN=filenamelen) :: outfile
  CHARACTER(LEN=filenamelen) :: stoke
  CHARACTER(LEN=filenamelen) :: dummy
  CHARACTER(LEN=filenamelen) :: store
  CHARACTER(LEN=filenamelen) :: filew,filea
  CHARACTER(LEN=filenamelen) :: noisefile
  CHARACTER(LEN=filenamelen) :: chline, chline1
  CHARACTER(LEN=*), PARAMETER :: code ="FastICA"
  CHARACTER(LEN=*), PARAMETER :: version = "1.0.0"

  type(paramfile_handle) :: handle

  REAL(kind=DP), DIMENSION(0:3) :: mono_dip
  REAL(kind=DP)                 :: theta_dip, phi_dip
  INTEGER(I4B)                  :: lowlreg

  !-------------------------------------------------------------------
  !                   get input parameter and create arrays
  !-------------------------------------------------------------------

  call date_and_time(values = values_time(:,1))
  if (nArguments() == 0) then
     paramfile=''
  else
     if (nArguments() /= 1) then
        print '("Usage: fastICA [parameter file name]")'
        stop 1
     endif
     call getArgument(1,paramfile)
  endif

  ! --- declaretion of intent
  PRINT*," "
  PRINT*,"               "//code//" "//version
  PRINT*," *** Component Separation with Fast ICA algorithm ***"
  PRINT*," "
  handle = parse_init(paramfile)

  ! --- first ask for total intensity or polarisation
  description = concatnl( &
       & " Do you want to analyse ",&
       & " 1) Total intensity ", &
       & " 2) Polarisation ")
  icatype = parse_int(handle, 'icatype', default=1, vmin=1, vmax=2, descr=description)
  if (icatype == 1) then 
     polarisation = .false.
     nsample = 1
     nstokes = 1
     qutype  = 0
  endif
  if (icatype == 2) then 
     polarisation = .true.
     description = concatnl( &
          & " Do you want to treat ", &
          & " 1) Q and U together (in the same vector)", &
          & " 2) Q and U separately (two different vectors)")
     qutype = parse_int(handle, 'qutype', default=1, vmin=1, vmax=2, descr=description)
     if (qutype == 1) then
        nsample = 2
        nstokes = 1
     endif
     if (qutype == 2) then
        nsample = 1
        nstokes = 2
     endif
  endif

  ! --- choose the number of input signals
  description = concatnl( &
       & " Enter the number of input signals (frequency channels): ")
  nsig = parse_int(handle, 'nsig', default=1, vmin=2, descr=description)
  
  ! --- gets parameters for ICA: max num of iteration, threshold for convergence,
  ! --- coefficients for Gauss and Tanh non-quadratic form, quadratic form selected
  description = concatnl( &
       & " Enter the maximum number of ICA iterations: ")
  maxNumIterations = parse_int(handle, 'n_iter', default = 20000, vmin = 1000, descr=description)
  description = concatnl( &
       & " Enter threshold for convergence")
  epsilon = parse_double(handle, 'epsilon', default = 1.0d-10, vmin=1.0d-10, descr=description)
  description = concatnl( &
       & " Enter a1 and a2 constants for Gauss and Tanh non-quadratic forms:")
  a1 = parse_double(handle, 'a1', default = 1.d0, vmin=1.d0, descr=description)
  a2 = parse_double(handle, 'a2', default = 1.d0, vmin=1.d0)
  description = concatnl( &
       & " Enter the kind of non-quadratic function", &
       & " p) power-law: resamble kurtosis", &
       & " g) gaussian", &
       & " t) Tanh")
  chline=''
  g = parse_string(handle, 'non_lin', default=chline, descr=description)
  g = TRIM(g)
  

  ! --- gets where you want to store data
  dummy = 'storedir'
  description = concatnl( &
       & " Enter the directory where you want to store data")
  store = parse_string(handle, dummy, default='.', descr=description)
  ! --- gets if noise realisations or rms maps are available
  description = concatnl( &
       & " Do you have noise realisation [true/false] ?")
  oknoise = parse_lgt(handle, 'oknoise', default=.false., descr=description)
  description = concatnl( &
       & " Do you have rms maps [true/false] ?")
  okrms = parse_lgt(handle, 'okrms', default=.false., descr=description)

  ! --- gets signal/noise/rms maps filename

  ALLOCATE(files(nstokes*nsample,nsig),noises(nstokes*nsample,nsig),rmsfiles(nstokes*nsample,nsig), stat=status)
  if (status /=0 ) call die_alloc(code,"files, noises, rmsfiles")
  chline=''
  ! --- now signals
  DO i = 1,nsig
     if (.not. polarisation) then
        write(output,"(i3)")i
        output=ADJUSTL(output)
        l=LEN_TRIM(output)
        dummy = 'file'//output(:l)
        description = concatnl( &
             & " Enter the # "//output(:l)//" signal file name (FITS file):")
        files(1,i) = parse_string(handle, dummy, default=chline, descr=description, filestatus='old')
        !
        ! --- noise files if needed
        !
        if (oknoise) then
           dummy = 'noise'//output(:l)
           description = concatnl( &
                & " Enter the # "//output(:l)//" noise file name (FITS file):")
           noises(1,i) = parse_string(handle, dummy, default=chline, descr=description, filestatus='old')
        endif
        !
        ! --- rms files if needed
        !
        if (okrms) then
           dummy = 'rms'//output(:l)
           description = concatnl( &
                & " Enter the # "//output(:l)//" rms file name (FITS file):")
           rmsfiles(1,i) = parse_string(handle, dummy, default=chline, descr=description, filestatus='old')
        endif
     else
        write(output,"(i3)")i
        output=ADJUSTL(output)
        l=LEN_TRIM(output)
        DO j = 1,nstokes*nsample
           write(stoke,"(i3)")j
           stoke=ADJUSTL(stoke)
           ll=LEN_TRIM(stoke)
           if (j == 1) dummy = 'fileQ'//output(:l)
           if (j == 2) dummy = 'fileU'//output(:l)
           description = concatnl( &
                & " Enter the "//stoke(:ll)//" Stoke parameter for # "//output(:l)//" signal file name (FITS file):")
           files(j,i) = parse_string(handle, dummy, default=chline, descr=description, filestatus='old')         
           !
           ! --- noise files if needed
           !
           if (oknoise) then
              if (j == 1) dummy = 'noiseQ'//output(:l)
              if (j == 2) dummy = 'noiseU'//output(:l)
              description = concatnl( &
                   & " Enter the "//stoke(:ll)//" Stoke parameter for # "//output(:l)//" noise file name (FITS file):")
              noises(1,i) = parse_string(handle, dummy, default=chline, descr=description, filestatus='old')
           endif
           !
           ! --- rms files if needed
           !
           if (okrms) then
              if (j == 1) dummy = 'rmsQ'//output(:l)
              if (j == 2) dummy = 'rmsU'//output(:l)
              description = concatnl( &
                   & " Enter the "//stoke(:ll)//" Stoke parameter for # "//output(:l)//" rms file name (FITS file):")
              rmsfiles(1,i) = parse_string(handle, dummy, default=chline, descr=description, filestatus='old')
           endif                     
        ENDDO
     endif
  enddo

  npixtot = getsize_fits(files(1,1), nmaps = nmaps, ordering = ordering, nside = nsmax, mlpol = mlpol)
  if (nsmax <= 0) then
     print*,"Keyword NSIDE not found in FITS header!"
     stop 1
  endif
  if (nsmax /= npix2nside(npixtot)) then
     print*,"FITS header keyword NSIDE does not correspond"
     print*,"to the size of the map!"
     stop 1
  endif

  ! --- check ordering scheme ---
  if ( (ordering /= 1) .and. (ordering /= 2)) then
     print*,"The ordering schem of the map must be RING or NESTED."
     print*,"No ordering specification is given in the FITS-header!"
     stop 1
  endif

  ! --- ask about monopole/dipole removal ---
  description = concatnl( &
       & "", &
       & " To improve the separation algorithm when working on the cut sky ", &
       & " you may want to remove the best fit monopole and dipole", &
       & " computed on the valid pixels out of the galactic cut.", &
       & " Do you want to :", &
       & " 0) Do the analysis on the raw map", &
       & " 1) Remove the monopole (mean) term ", &
       & " 2) Remove the monopole and dipole terms ")
  lowlreg = parse_int(handle, 'regression', vmin=0, vmax=2, default=0, descr=description)

  !-------------------------------------------------
  !       ask for noise in the data
  !------------------------------------------------
  description = concatnl( &
       & "", &
       & " Noise in the data [true/false]? ")
  noise_in_data = parse_lgt(handle, 'noise_in_data', default=.false., descr=description)
  if (noise_in_data) then
     description = concatnl( &
          & "", &
          & " Do you want to: ", &
          & " 1) Account for noise in data", &
          & " 2) Leave data as they are")
     cure_noise = parse_int(handle, 'cure_noise', vmin = 1, vmax = 2, default = 1, descr=description)
     if (cure_noise == 1) then 
        description = concatnl( &
             & "", &
             & " Do you want to: ", &
             & " 1) Use Nobs map plus noise per observation", &
             & " 2) Use a constant guess rms ")
        choice_rmsprior = parse_int(handle, 'rmsprior', vmin = 1, vmax = 2, default = 1, descr=description)
     endif
     if (choice_rmsprior == 1) then

        ALLOCATE(rms_guess(nstokes*nsample,nsig),stat=status)
        if (status /= 0) call die_alloc(code,"rms_guess")

        do i = 1, nsig
           write(output,"(i3)")i
           output=ADJUSTL(output)
           l=LEN_TRIM(output)
           if (.not. polarisation) then
              dummy = "rms"//output(:l)
              description = concatnl( &
                   & "", &
                   & " Enter the rms per observation for channel # "//output(:l))
              rms_guess(1,i) = parse_double(handle, dummy, vmin= 0.0_dp, default = 0.0_dp, descr=description)
           else
              do j=1,nstokes*nsample
                 write(stoke,"(i3)")j
                 stoke=ADJUSTL(stoke)
                 ll=LEN_TRIM(stoke)
                 if (j == 1) dummy = "rmsQ"//output(:l)
                 if (j == 2) dummy = "rmsU"//output(:l)
                 description = concatnl( &
                      & "", &
                      & " Enter the rms per observation for the "//stoke(:ll)//" Stokes parameter for channel # "//output(:l))
                 rms_guess(j,i) = parse_double(handle, dummy, vmin= 0.0_dp, default = 0.0_dp, descr=description)
              enddo
           endif
        enddo
     endif
     if (choice_rmsprior == 2) then
        do i = 1, nsig
           write(output,"(i3)")i
           output=ADJUSTL(output)
           l=LEN_TRIM(output)
           if (.not. polarisation) then
              dummy = "rms"//output(:l)
              description = concatnl( &
                   & "", &
                   & " Enter the guess rms for channel # "//output(:l))
              rms_guess(1,i) = parse_double(handle, dummy, vmin = 0.0_dp, default = 0.0_dp, descr=description)
           else
              do j=1,nstokes*nsample
                 write(stoke,"(i3)")j
                 stoke=ADJUSTL(stoke)
                 ll=LEN_TRIM(stoke)
                 if (j==1) dummy="rmsQ"//output(:l)
                 if (j==2) dummy="rmsU"//output(:l)
                 description = concatnl( &
                      & "", &
                      & " Enter the guess rms for the "//stoke(:ll)//" Stokes parameters for the channel # "//output(:l))
                 rms_guess(j,i) = parse_double(handle,dummy, vmin= 0.0_dp, default = 0.0_dp, descr=description)
              enddo
           endif
        enddo
     endif
     
  endif

  !-------------------------------------------------------------
  !            ask for sky cut and read the input signals
  !-------------------------------------------------------------
  ALLOCATE(mask(0:npixtot-1,1:nmaps),stat=status)
  if (status /= 0) call die_alloc(code,"mask")
  description = concatnl( &
       & "", &
       & " Do you want :", &
       & " 1) Full-sky analysis", &
       & " 2) |b| > b_cut", &
       & " 3) b_min < b < b_max", &
       & " 4) b_min < |b| < b_max with b_min and b_max >= 0", &
       & " 5) custom cut")
  choice_cut = parse_int(handle, 'cut', vmin=1, vmax=5, default=1, descr=description)

  ! --- a real cut is here: symmetric in general
  if (choice_cut /= 1) PRINT *,"      "//code//"> Computing/Inputting Mask file"
  if (choice_cut >=2 .and. choice_cut <= 4) then
     description = concatnl( &
          & " Enter b_min and b_max [-90,90]:")
     b_min = parse_double(handle, 'bmin', vmin=-90.0_dp, vmax=90.0_dp, default=0.0_dp, descr=description)
     b_max = parse_double(handle, 'bmax', vmin=-90.0_dp, vmax=90.0_dp, default=0.0_dp)
     description = concatnl( &
          & " Enter l_min and l_max [0,360]:")
     l_min = parse_double(handle, 'lmin', vmin=0.0_dp, vmax=360.0_dp, default=0.0_dp, descr=description)
     l_max = parse_double(handle, 'lmax', vmin=0.0_dp, vmax=360.0_dp, default=0.0_dp, descr=description)
     call get_gen_cut(nmaps, ordering, npixtot, choice_cut, l_min, l_max, b_min, b_max, mask, ncount)
  endif

  ! --- full-sky analysis: cut is zero
  if (choice_cut == 1) then
     mask(:,1:nmaps) = 1
     ncount = npixtot
  else if (choice_cut == 5) then
     ! --- custom cut: user supplied FITS file
     ! --- first the mask file
     description = concatnl( &
          & " Enter mask file name (FITS file): ")
     maskfile = parse_string(handle, 'maskfile', default=chline, descr=description, filestatus='old')
     call input_map(maskfile, mask(0:,1:nmaps), npixtot, nmaps, fmissval=fmissval)     
  endif

  if (choice_cut /= 1) then
     ! --- count pixels after galactic cut 
     ncount = 0
     ncut = 0
     do i = 0,npixtot-1
        if (mask(i,1) == 1) then
           ncount = ncount + 1
        else
           ncut = ncut + 1
        endif
     enddo
     ! --- extract index for non-zero pixel in mask
     !
     ALLOCATE(left_pixel(1:ncount),stat=status)
     if (status /= 0) call die_alloc(code,"left_pixel")
     
     ALLOCATE(cut_pixel(1:ncut),stat=status)
     if (status /= 0) call die_alloc(code,"cut_pixel")
     
     k = 1
     kk = 1
     do j = 1, nmaps
        do i = 1,npixtot
           if (mask(i-1,j) == 1.) then
              left_pixel(k) = i-1
              k = k + 1
           else
              cut_pixel(kk) = i-1
              kk = kk + 1
           endif
        enddo
     enddo
     PRINT*,"      "//code//"> Total number of pixel after cut  = ",k-1
     PRINT*,"      "//code//"> Total number of pixel in the cut = ",kk-1 
  endif

  ! --- now enlarge dimension for arrays in case of polarisation and 
  ! --- separation performed with Q and U together
  if (polarisation) ncount = ncount * nsample

  fmissval = 0.0_sp

!-------------------------------------------------------------
!                allocate space for arrays
!-------------------------------------------------------------

  ALLOCATE(mapIO(0:npixtot-1,1),err(0:npixtot-1,1),stat=status)
  if (status /= 0) call die_alloc(code,"mapIO or err")

  ALLOCATE(sig(nsig,ncount),xsig(nsig,ncount),n(nsig,ncount),errmap(nsig,ncount),stat=status)
  if (status /= 0) call die_alloc(code,"sig, xsig, n, errmap")

  ALLOCATE(med(nsig), wmat(nsig,nsig), dwmat(nsig,nsig), stat=status)
  if (status /= 0) call die_alloc(code,"med, wmat, dwmat, rms")

  ALLOCATE(s(nsig,nsig),a(nsig,nsig),w(nsig,nsig),stat=status)
  if (status /= 0) call die_alloc(code,"s, a, w")

  !----------------------------------------------------------
  !            input signal maps
  !----------------------------------------------------------
  PRINT*,"      "//code//"> Inputting signal maps"

  do j = 1,nstokes
     do k = 1,nsample
        if (qutype == 1) then
           PRINT*,"      "//code//"> Working on the "//char(48+k)//"-th sample"
        else if (qutype == 2) then
           PRINT*,"      "//code//"> Working on the "//char(48+j)//"-th Stokes parameter"
        endif
        do i = 1,nsig
           call input_map(files(j*k,i), mapIO(0:,1:1), npixtot, nmaps, fmissval=fmissval) 
!           print*,mapIO(0:9,1:1)
!           print*," "
           !----------------------------------------------------------------
           !                 remove dipole
           !----------------------------------------------------------------
           
           if (lowlreg > 0) then
              PRINT *,"      "//code//"> Remove monopole (and dipole) from input maps"
              if (choice_cut == 1) then
                 cos_theta_cut = cos(1.0_dp)
                 zbounds = (/cos_theta_cut, -cos_theta_cut/)
!                 call remove_dipole(nsmax, mapIO(0:npixtot-1,1), ordering, lowlreg, &
!                      & mono_dip(0:), zbounds, fmissval=fmissval,mask=mask)
              else
!                 call remove_dipole_custom(nsmax, mapIO(0:npixtot-1,1), 1, ordering, lowlreg, &
!                      & mono_dip(0:), cut_pixel, fmissval)
              endif
              
              if (fmissval /= 0.0) then
                 do m = 1,nmaps
                    do ll = 0,npixtot-1 
                       if (abs(mapIO(ll,m)/fmissval-1.0) < 1.e-6*fmissval) then
                          mapIO(ll,m) = 0.0_sp
                       endif
                    enddo
                 enddo
              endif
              
              write(unit=*,fmt="(a,g13.3)")  " Monopole = ",mono_dip(0)
              if (lowlreg > 1) then
                 call vec2ang( mono_dip(1:3), theta_dip, phi_dip)
                 write(unit=*,fmt="(a,g13.3)")  " Dipole   = ", &
                      & sqrt(sum(mono_dip(1:3)**2))
                 write(unit=*,fmt="(a,g10.3,', ',g10.3,a)") &
                      & "(long.,lat.) = (",phi_dip/PI*180.0,90.0-theta_dip/PI*180.0,") Deg"
              endif
              
        
           endif  ! --- done with removing dipole

           do m = 1,ncount/nsample
              if (choice_cut == 1) then
                 sig(i,m+(k-1)*ncount/nsample) = mapIO(m-1,1)
              else
                 sig(i,m+(k-1)*ncount/nsample) = mapIO(left_pixel(m),1) 
              endif
           enddo

        enddo
     enddo
     ! put noise matrix to zero if noise is not considered
     do i=1,nsig
        do m=1,nsig
           s(i,m) = 0.d0
        enddo
     enddo

     !----------------------------------------------------------
     !           input Nobs maps
     !----------------------------------------------------------
     IF (cure_noise == 1 .and. choice_rmsprior == 1) then
        PRINT*,"      "//code//"> Inputting Nobs maps"
        do i = 1,nsig
           s = 0.0_dp
           call input_map(noises(j*k,i), mapIO(0:,1:nmaps), npixtot, nmaps, fmissval = fmissval)
           do m = 1,ncount/nsample
              if (choice_cut == 1) then
                 n(i,m+(k-1)*ncount/nsample) = mapIO(m-1,1)
              else
                 n(i,m+(k-1)*ncount/nsample) = mapIO(left_pixel(m),1)
              endif
           enddo
           med(i) = SUM(n(i,:))/ncount
           ! ------------------------------------------------
           ! Attention!!!! Depending on the S/N ratio of your
           ! data you would ask for an rms "uniformied"
           !-------------------------------------------------
           ! --- in this way the rms is "uniformied"
           rms_guess(j*k,i) = rms_guess(j*k,i)/sqrt(med(i))
           s(i,i) = rms_guess(j*k,i)**2
        enddo
        ! --- put the non-diagonal component of s equal to 0
        do i = 1,nsig
           do m = 1,nsig
              if (i /= m) s(i,m) = 0.d0
           enddo
           PRINT*,"      "//code//"Diagonal of the noise matrix:"
           PRINT*,s(i,i)
        enddo
     endif
     
     if (cure_noise == 1 .and. choice_rmsprior == 2) then
        do i = 1,nsig
           do m = 1,nsig
              if (i == m) s(i,i) = rms_guess(j*k,i)**2
              if (i /= m) s(i,m) = 0.d0
           enddo
           PRINT*,"      "//code//"Diagonal of the noise matrix:"
           PRINT*,s(i,i)
        enddo
     endif
     
     
     if (lowlreg > 0) then
        header = ''
        PRINT*,"      "//code//"> Writing input map with modes removed to FITS files"
        call add_card(header,'COMMENT','------------------------------------------------')
        call add_card(header,'COMMENT','      Sky Map Pixelisation Specific Keywords    ')
        call add_card(header,'COMMENT','------------------------------------------------')
        call add_card(header,'PIXTYPE','HEALPIX','HEALPIX Pixelisation')
        if (ordering == 1) call add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED ')
        if (ordering == 2) call add_card(header,'ORDERING','NESTED',  'Pixel ordering scheme, either RING or NESTED ')
        call add_card(header,'NSIDE'   ,nsmax,  'Resolution parameter for HEALPIX')
        call add_card(header,'FIRSPIX',0,'First pixel # (0 based)')
        call add_card(header,'LASTPIX',npixtot-1,'Last pixel # (0 based)')
        call add_card(header)
        
        call add_card(header,'COMMENT','------------------------------------------------')
        call add_card(header,'COMMENT','      Planck Simulation Specific Keywords       ')
        call add_card(header,'COMMENT','------------------------------------------------')
        call add_card(header,'EXTNAME','ICADATA')
        call add_card(header,'CREATOR',code,        'Software creating the FITS file')
        call add_card(header,'VERSION',version,     'Version of the analysis software')
        call add_card(header,'COMMENT','------------------------------------------------')
        call add_card(header,'COMMENT','      Data Destription Specific Keywords       ')
        call add_card(header,'COMMENT','------------------------------------------------')
        call add_card(header,'PDMTYPE',  'COMPMAP','Planck Data Model Type')
        
        nlheader = SIZE(header) 
        do i = 1,nsig
           mapIO(:,1) = 0.
           do k = 1,nsample
              do m = 1,ncount/nsample
                 if (choice_cut == 1) then
                    mapIO(m-1,1) = sig(i,m+(k-1)*ncount/nsample)
                 else
                    mapIO(left_pixel(m-1),1) = sig(i,m+(k-1)*ncount/nsample)
                 endif
              enddo
              where (mapIO(:,1) == 0.)
                 mapIO(:,1) = -1.6375e30
              end where
              write(output,"(i3)")i
              output=ADJUSTL(output)
              l=LEN_TRIM(output)
              if (.not. polarisation) then
                 outfile = trim(store)//"/x"//output(:l)//".fits"
              else
                 if (j == 1 .or. k == 1) outfile = trim(store)//"/xQ"//output(:l)//".fits"
                 if (j == 2 .or. k == 2) outfile = trim(store)//"/xU"//output(:l)//".fits"
              endif
              call write_bintab(mapIO, npixtot, nmaps ,header, nlheader, outfile)
           enddo
        enddo
     endif
     
     ! --- now subtract the mean from input signals ---
     do i = 1,nsig
        med(i) = SUM(sig(i,:))/ncount
        sig(i,:) = sig(i,:) - med(i)
     enddo
     
     ! --- these are the actual data on which the algorithm works
     xsig = sig
     
     PRINT*,"      "//code//"> Begin whitening "
     call whitening(code, nsig, ncount, xsig, s, wmat, dwmat) 
     
     header = ''
     PRINT*,"      "//code//"> Writing whitened map to FITS files"
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','      Sky Map Pixelisation Specific Keywords    ')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'PIXTYPE','HEALPIX','HEALPIX Pixelisation')
     if (ordering == 1) call add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED ')
     if (ordering == 2) call add_card(header,'ORDERING','NESTED',  'Pixel ordering scheme, either RING or NESTED ')
     call add_card(header,'NSIDE'   ,nsmax,  'Resolution parameter for HEALPIX')
     call add_card(header,'FIRSPIX',0,'First pixel # (0 based)')
     call add_card(header,'LASTPIX',npixtot-1,'Last pixel # (0 based)')
     call add_card(header)
     
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','      Planck Simulation Specific Keywords       ')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'EXTNAME','ICADATA')
     call add_card(header,'CREATOR',code,        'Software creating the FITS file')
     call add_card(header,'VERSION',version,     'Version of the analysis software')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','      Data Destription Specific Keywords       ')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','Whitened map')
     call add_card(header,'PDMTYPE',  'COMPMAP','Planck Data Model Type')
     
     nlheader = SIZE(header) 
     do i = 1,nsig
        mapIO(:,1) = 0.
        do k=1,nsample
           do m = 1,ncount/nsample
              if (choice_cut == 1) then
                 mapIO(m-1,1) = xsig(i,m+(k-1)*ncount/nsample)
              else
                 mapIO(left_pixel(m),1) = xsig(i,m+(k-1)*ncount/nsample)
              endif
           enddo
           write(output,"(i3)")i
           output=ADJUSTL(output)
           l=LEN_TRIM(output)
           if (.not. polarisation) then
              outfile = trim(store)//"xw"//output(:l)//".fits"
           else
              if (j == 1 .or. k == 1) outfile = trim(store)//"xwQ"//output(:l)//".fits"
              if (j == 2 .or. k == 2) outfile = trim(store)//"xwU"//output(:l)//".fits"
           endif
!           write(*,*) 'sono qui con ',trim(outfile)
!           call write_bintab(mapIO, npixtot, nmaps, header, nlheader,outfile)
!           write(*,*) ' fatto'
        enddo
     enddo
     ! ------------------------------------------------------
     !                 begin with ICA
     ! ------------------------------------------------------
     PRINT*,"      "//code//"> Begin FastICA "
     call fpica(code, nsig, ncount, xsig, s, wmat, dwmat, g, a1, a2, epsilon, maxNumIterations, a, w, conv)
     
     if (noise_in_data .or. oknoise) then
        PRINT*,"      "//code//"> Trying to estimate noise"
        errmap=MATMUL(w,n)
     endif
     
     !-----------------------------------------------------------------------
     !     output is linear combination of inputs acccording to w matrix
     !-----------------------------------------------------------------------
     xsig = MATMUL(w,sig)
     PRINT*,"      "//code//"> Normalize outputs and noise to the "//char(48+1)//"-th channel"
     
     do i = 1,nsig
!        xsig(i,:) = xsig(i,:) * a(1,i) 
        errmap(i,:) = errmap(i,:) * a(1,i) 
     end do
     
     PRINT*,"      "//code//"> Save results on files "
     
     !-------------------------------------------------------------
     !         write the reconstructed components to FITS files
     !-------------------------------------------------------------
     header = ''
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','      Sky Map Pixelisation Specific Keywords    ')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'PIXTYPE','HEALPIX','HEALPIX Pixelisation')
     if (ordering == 1) call add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED ')
     if (ordering == 2) call add_card(header,'ORDERING','NESTED',  'Pixel ordering scheme, either RING or NESTED ')
     call add_card(header,'NSIDE'   ,nsmax,  'Resolution parameter for HEALPIX')
     call add_card(header,'FIRSPIX',0,'First pixel # (0 based)')
     call add_card(header,'LASTPIX',npixtot-1,'Last pixel # (0 based)')
     call add_card(header)
     
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','      Planck Simulation Specific Keywords       ')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'EXTNAME','ICADATA')
     call add_card(header,'CREATOR',code,        'Software creating the FITS file')
     call add_card(header,'VERSION',version,     'Version of the analysis software')
     call add_card(header,'PDMTYPE',  'COMPMAP','Planck Data Model Type')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'COMMENT','      FastICA Specific Keywords                 ')
     call add_card(header,'COMMENT','------------------------------------------------')
     call add_card(header,'NONLINE',TRIM(g),           'Non-quad. form used by FastICA')
     if (choice_cut >= 2 .and. choice_cut <= 4) then
        call add_card(header,'BMINCUT',b_min,           'b_min for cut')
        call add_card(header,'BMAXCUT',b_max,           'b_max for cut')
        call add_card(header,'LMINCUT',l_min,           'l_min for cut')
        call add_card(header,'LMAXCUT',l_max,           'l_max for cut')
     endif
     if (choice_cut == 1) then
        call add_card(header,'COMMENT','Full-sky analysis')
     endif
     if (choice_cut == 5) then
        call add_card(header,'COMMENT','Mask file = '//TRIM(maskfile))
     endif
     call add_card(header,'CONVERGE',conv,            '0=non-conv,1=conv') 
     
     nlheader = SIZE(header)
     do i = 1,nsig
        do k = 1,nsample
           do m = 1,ncount/nsample
              if (choice_cut == 1) then
                 mapIO(m-1,1) = xsig(i,m+(k-1)*ncount/nsample)
                 if (oknoise .or. cure_noise == 1) err(m-1,1) = errmap(i,m+(k-1)*ncount/nsample)
              else
                 mapIO(left_pixel(m),1) = xsig(i,m+(k-1)*ncount/nsample)
                 if (oknoise .or. cure_noise == 1) err(left_pixel(m),1) = errmap(i,m+(k-1)*ncount/nsample)
              endif
           enddo
	do ipix=0,npixtot-1 
           if (mapIO(ipix,1) == 0.) mapIO(ipix,1)=-1.6375e+30
        enddo
        print*,' ho finito'
           write(output,"(i3)")i
           output = ADJUSTL(output)
           l=LEN_TRIM(output)
           if (.not. polarisation) then
              outfile = trim(store)//"/u"//output(:l)//".fits"
              noisefile = trim(store)//"/n"//output(:l)//".fits"
           else
              if (j==1 .or. k == 1) then
                 outfile = trim(store)//"/u"//output(:l)//"_Q.fits"
                 noisefile = trim(store)//"/n"//output(:l)//"_Q.fits"
              endif
              if (j==2 .or. k == 2) then
                 outfile = trim(store)//"/u"//output(:l)//"_U.fits"
                 noisefile = trim(store)//"/n"//output(:l)//"_U.fits"
              endif
           endif
           call write_bintab(mapIO, npixtot, nmaps, header, nlheader, outfile)
           if (oknoise .or. cure_noise == 1) then 
              call write_bintab(err, npixtot, nmaps, header, nlheader, noisefile)
           endif
        enddo
     enddo

     !----------------------------------
     ! write files with W and A matrices
     !----------------------------------
     if (.not. polarisation) then
        filew = trim(store)//"/w.txt"
        filea = trim(store)//"/a.txt"
     else
        do k = 1,nsample
           if (j == 1 .or. k == 1) then
              filew = trim(store)//"/w_Q.txt"
              filea = trim(store)//"/a_Q.txt"
           endif
           if (j == 2 .or. k == 2) then
              filew = trim(store)//"/w_U.txt"
              filea = trim(store)//"/a_U.txt"
           endif
        enddo
     endif
     OPEN(1,file=filew,status="replace")
     OPEN(2,file=filea,status="replace")
     do i=1,nsig
        write(1,*) (w(i,m),m=1,nsig)
        write(2,*) (a(i,m),m=1,nsig)
     enddo
     CLOSE(2)
     CLOSE(1)
  enddo
  !-------------------------------------------------
  !           deallocate memory for arrays
  !-------------------------------------------------
  deallocate(sig, xsig, med, wmat, dwmat, s, a, w, mapIO)
  deallocate(mask, err, errmap)
  if (choice_cut /= 1) deallocate(left_pixel, cut_pixel)
  
  !--------------------------------------
  !          report card
  !--------------------------------------
  call date_and_time(values = values_time(:,2))
  
  values_time(:,1) = values_time(:,2) - values_time(:,1)
  clock_time =  (  (values_time(3,1)*24 &
       &           + values_time(5,1))*60. &
       &           + values_time(6,1))*60. &
       &           + values_time(7,1) &
       &           + values_time(8,1)/1000.

  WRITE(*,9000) " "
  WRITE(*,9000) " Report Card for "//code//" analysis run"
  WRITE(*,9000) "----------------------------------------"
  WRITE(*,9000) " "
  WRITE(*,9000) "       Temperature alone"
  WRITE(*,9000) " "
  WRITE(*,9010) " Number of Signals            : ",nsig
  WRITE(*,9010) " Number of pixels             : ",npixtot
  WRITE(*,9010) " Max Num of Iteration for ICA : ",maxNumIterations
  WRITE(*,9020) " Threshold for convergence    : ",epsilon
  WRITE(*,9000) " Non-quadratic form adopted   : "//TRIM(g)
  if (g .ne. 'p') then
     WRITE(*,9020) " a1                          : ",a1
     WRITE(*,9020) " a2                          : ",a2
  endif
  WRITE(*,9000)"       Input files"
  do i=1,nsig
     if (.not. polarisation) then
        WRITE(*,9001) " signal",i,"                     : "//TRIM(files(1,i)) 
        if (oknoise) then
           WRITE(*,9001)" noise",i,"                     :"//TRIM(noises(1,i))
        endif
        if (okrms) then
           WRITE(*,9001)" rms",i,"                     :"//TRIM(rmsfiles(1,i))
        endif
     else
        WRITE(*,9001) "signal Q",i,"                     : "//TRIM(files(1,i))
        WRITE(*,9001) "signal U",i,"                     : "//TRIM(files(2,i))
        if (oknoise) then
           WRITE(*,9001) "noise Q",i,"                     : "//TRIM(noises(1,i))
           WRITE(*,9001) "noise U",i,"                     : "//TRIM(noises(2,i))
        endif
        if (okrms) then
           WRITE(*,9001) "rms Q",i,"                     : "//TRIM(rmsfiles(1,i))
           WRITE(*,9001) "rms U",i,"                     : "//TRIM(rmsfiles(2,i))
        endif
     endif
  enddo
  WRITE(*,9020) " Total clock time [s] : ",clock_time

  !-----------------------------------
  !           end of routine
  !-----------------------------------
  PRINT*,"         "//code//"> Normal completion "

9000 format(a)
9001 format(a,i2,a)
9010 format(a,i16)
9020 format(a,g20.5)

  STOP
END PROGRAM fastica






