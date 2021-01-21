module co_util

    implicit none

contains
    recursive subroutine co_bcast(coarray, source_image, first_image, last_image)
        implicit none
#if NO_COARRAYS
        real(kind=8), intent(inout) :: coarray(:,:,:,:)
#else
        real(kind=8), intent(inout) :: coarray(:,:,:,:)[*]
#endif
        integer, intent(in) :: source_image, first_image, last_image
        integer :: dest_image
        integer :: x,y
#ifndef NO_COARRAYS
        if (first_image==last_image) return

        if (source_image/=first_image) then
            dest_image=first_image

            if (this_image()==source_image) then
                do y=1,size(coarray,4)
                    do x=1,size(coarray,3)
                        coarray(:,:,x,y)[dest_image] = coarray(:,:,x,y)
                    enddo
                enddo
                sync images(dest_image)
            elseif (this_image()== dest_image) then
                sync images(source_image)
            endif

            if (this_image()<source_image) then
                call co_bcast(coarray, dest_image, dest_image, source_image-1)
            else
                call co_bcast(coarray, source_image, source_image, last_image)
            endif
        else
            dest_image = ((last_image-first_image)+1)/2 + first_image

            if (this_image()==source_image) then
                do y=1,size(coarray,4)
                    do x=1,size(coarray,3)
                        coarray(:,:,x,y)[dest_image] = coarray(:,:,x,y)
                    enddo
                enddo
                sync images(dest_image)
            elseif (this_image()== dest_image) then
                sync images(source_image)
            endif

            if (this_image()<dest_image) then
                call co_bcast(coarray, source_image, first_image, dest_image-1)
            else
                call co_bcast(coarray, dest_image, dest_image, last_image)
            endif
        endif
#endif
    end subroutine co_bcast

end module co_util

 ! in case both images have to hit the exact same sync point
! if ((this_image()==source_image).or.(this_image()==dest_image)) then
!     if (this_image()==source_image) then
!         sync_image = dest_image
!     else
!         sync_image = source_image
!     endif
!     sync images(sync_image)
! endif
