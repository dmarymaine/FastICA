module utils
  ! Author: Davide Maino
  ! CHANGES:
  ! original version 27/FEB/2001 by D.Maino @MPA
  ! modified 16/MAR/2001: add get_mask_cut 
  ! modified 16/JAN/2004: leave only get_gen_cut for HEALPix
  
  use healpix_types
  use pix_tools, ONLY : pix2ang_ring, pix2ang_nest

contains

  !==========================================================================
  subroutine get_gen_cut(nmaps, ordering, npixtot, choice_cut, l_min, l_max, b_min, b_max, mask, ncount)
    !========================================================================
    implicit none
    integer(i4b) :: ncount,nside,npixtot, i,choice_cut,ordering, nmaps
    real(sp), dimension(1:npixtot,1:nmaps) :: mask
!    real(sp), dimension(0:2) :: vec
    real(dp) :: l_min,l_max,b_min,b_max,l_mi,l_ma,b_mi,b_ma,theta,phi
!,norm_inv

    nside = sqrt(npixtot/12.0_sp)

    l_mi = l_min * pi /180.0_dp
    l_ma = l_max * pi /180.0_dp

    if (b_min >= 0.0) b_mi = halfpi - b_min * pi/180.0_dp
    if (b_min < 0.0)  b_mi = halfpi + abs(b_min) * pi/180.0_dp 
    if (b_max >= 0.0) b_ma = halfpi - b_max * pi/180.0_dp
    if (b_max < 0.0)  b_ma = halfpi + abs(b_max) * pi/180.0_dp

    ncount = 0 
    DO i=0,npixtot-1
       if (ordering == 1) then
          call pix2ang_ring(nside,i,theta,phi)
       else 
          call pix2ang_nest(nside,i,theta,phi)
       endif
   
       if (phi>=l_mi .and. phi <=l_ma) then
          if(choice_cut == 2 ) then
             if(theta .le. b_ma .or. theta .ge. b_mi) then
                mask(i+1,1:nmaps) = 1.
                ncount = ncount + 1
             endif
          else if(choice_cut == 3) then
             if(theta .ge. b_ma .and. theta .le. b_mi) then
                mask(i+1,1:nmaps) = 1.
                ncount = ncount + 1
             endif
          else if (choice_cut == 4) then
             if( (theta.ge.b_ma .and. theta.le.b_mi).or.(theta.ge.(b_ma+halfpi).and. theta.le.(b_mi+halfpi))) then
                mask(i+1,1:nmaps) = 1.
                ncount = ncount +1
             endif
          endif
       endif
    ENDDO
    RETURN
  end subroutine get_gen_cut

end module utils
