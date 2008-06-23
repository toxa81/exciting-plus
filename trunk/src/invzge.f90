      subroutine invzge(A,ndim)                                                                                                        
      implicit none                                                                                                                               
! passed var                                                                                                                                      
      integer ndim                                                                                                                                
      complex*16 A(ndim,ndim)                                                                                                                     
! local var                                                                                                                                       
      integer lwork,nb,info                                                                                                                       
      real*8 ,allocatable :: work(:)                                                                                                              
      integer ,allocatable :: ipiv(:)                                                                                                             
      integer i,j                                                                                                                                 
      integer, external :: ilaenv                                                                                                                 
      
      
      nb = ilaenv(1,'zgetri','U',ndim,-1,-1,-1)                                                                                                   
      lwork = ndim * nb                                                                                                                           
      allocate(work(2*lwork),ipiv(ndim))                                                                                                          
                                                                                                                                                  
      call zgetrf(ndim,ndim,A,ndim,ipiv,info)                                                                                           
      if(info.ne.0) then                                                                                                                            
        write(*,*)'inverse_he_matrix: error factorization'
	stop
      endif                                                                                 
                                                                                                                                                  
      call zgetri(ndim,A,ndim,ipiv,work,lwork,info)                                                                                                 
      if(info.ne.0) then                                                                                                                            
        write(*,*)'inverse_he_matrix: error inversion'
	stop
      endif                                                                                     
                                                                                                                                                  
      deallocate(work,ipiv)                                                                                                                       
                                                                                                                                                  
      end                  
