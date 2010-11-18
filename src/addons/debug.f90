subroutine dbg_open_file
use modmain
implicit none
open(fdbgout,file=trim(adjustl(fdbgname)),form="FORMATTED",status="OLD",&
  position="APPEND")
return
end

subroutine dbg_close_file
use modmain
implicit none
close(fdbgout)
return
end
