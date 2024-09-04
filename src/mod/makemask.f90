! -*- F90 -*-
program makemask
  use healpix_types
  use fitstools
  use head_fits

  implicit none
  integer(i4b) :: ipix,nside,npixtot,nlheader,npix,status
  integer(i4b), dimension(:), allocatable :: ind

  real(sp), dimension(:,:), allocatable :: mask
  real(sp), PARAMETER :: blank = -1.6375e30

  logical(lgt) :: there

  character(len=filenamelen) :: infile,maskfile
  character(len=80), dimension(1:120) :: header
  character(len=2) :: card
  character(len=8), PARAMETER :: code = 'makeMask'

! -*- Begin
1 print*,'  '//code//'> enter the pixel file name:'
  read(*,'(a)') infile
  inquire(file=infile,exist=there)
  if (.not.there) then
     print*,' '//code//'> input file not found!'
     goto 1
  endif

  print*,' '//code//'> enter the total number of pixels in the file'
  read(*,*) npix
  allocate(ind(npix),stat=status)
  if (status /= 0) then
     stop 'cannot allocate ind'
  endif
  
  print*,' '//code//'> enter the nside of the corresponding pixels numbers'
  read(*,*) nside
  npixtot = 12*nside**2
  write(card,'(i2)') nside
  maskfile = 'mask_h'//trim(card)//'_ring.fits'
  print*,' '//code//'> output file is ',maskfile

  print*,' '//code//'> Reading pixel file...'
  open(1,file=infile,status='old')
  read(1,*) ind
  close(1)

  allocate(mask(0:npixtot-1,1),stat=status)
  if(status /= 0) then
     stop 'cannot allocate mask'
  endif
  mask = 0
  
  mask(ind(:),1) = 1
  where (mask(:,1) == 0) 
     mask(:,1) = blank
  end where
  deallocate(ind)

  header = ''
  print*,' '//code//'> Writing mask to FITS file'
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'COMMENT','     Sky Map Pixelisation Specific Keywords    ')
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'PIXTYPE','HEALPIX','HEALPIX pixelisation')
  call add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED')
  call add_card(header,'NSIDE'   ,nside,   'Resolution parameter for HEALPIX')
  call add_card(header,'FIRSTPIX',0,'First pixel # (0 based)')
  call add_card(header,'LASTPIX',npixtot-1,'Last pixel # (0 based)')
  call add_card(header)

  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'COMMENT','     Planck Simulation Specific Keywords       ')
  call add_card(header,'COMMENT','-----------------------------------------------')
  call add_card(header,'EXTNAME','SIMULATION')
  call add_card(header,'CREATOR',code,       'Software creating the FITS file')
  call add_card(header,'DMRDMTY','MASK',      'DMR data model type')

  nlheader=size(header)
  call write_bintab(mask(0,1),npixtot,1,header,nlheader,maskfile)

  deallocate(mask)
end program makemask


