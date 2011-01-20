module mod_madness

integer m_ngvec
integer, allocatable :: m_ngknr(:)
integer, allocatable :: m_igkignr(:,:)

complex(8), allocatable :: m_wann_unkmt(:,:,:,:,:)
complex(8), allocatable :: m_wann_unkit(:,:,:)

real(8) m_wanpos(3)

end module