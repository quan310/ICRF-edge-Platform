Module FileIO
    Use TypeModule
    Use ParticleModule  
    implicit none
    contains   
    subroutine LoadParticle(InputSpecy,InputParticleBundle,Status)
        implicit none
        Type(OneSpecy),intent(in) ::  InputSpecy
        Type(ParticleBundle),intent(out) :: InputParticleBundle
        Logical,intent(out) ::  Status
        Logical :: alive
        Character(len=20) :: Filename
        Integer(4) :: i
        
        Write(filename,*) trim(InputSpecy%Name),".dat"
        Inquire(file=filename,exist=alive)
        If(alive) then
               Open (10,file=filename)
               InputParticleBundle%Name=InputSpecy%Name
               Read(10,*) InputParticleBundle%NPar
               Write(*,*) 'Loading ', trim(InputParticleBundle%Name), InputParticleBundle%NPar,' Please Wait...' 
               Write(filename,*) trim(InputSpecy%Name),".dat"
               do i=1,InputParticleBundle%NPar
                     Read(10,*)  InputParticleBundle%PhaseSpace(i)
               end do
               Status=.False.
               Write(*,*) 'Loading ', trim(InputParticleBundle%Name),' Complete!'   
         else
               Write(*,*) 'Can not find the file for ', trim(InputSpecy%Name),' the particle will be randomly initilalized.'    
               Status=.True. 
         End if        
        return
    end subroutine  LoadParticle      
    subroutine DumpGrid(InputGrid)
        Implicit none
        Type(Grid),intent(in) :: InputGrid
        Character(len=20) :: Filename
        Integer(4):: i,j,Index
        write(filename,*) trim(InputGrid%Name),".dat"
        Write(*,*) "Saving ",trim(InputGrid%Name)," Please wait..."
        open (10,file=filename)
        Write(10,*)  InputGrid%Nx,InputGrid%Ny
        Write(10,*)  InputGrid%Dx,InputGrid%Dy
        Index=0
        do i=1,InputGrid%Ny
               Write(10,FMt="(<InputGrid%Nx>D23.15)")  (InputGrid%Value(Index+j),j=1,InputGrid%Nx)
               Index=Index+InputGrid%Nx 
        end do
        close(10)
        Write(*,*) "Save ",trim(InputGrid%Name),"Complete!"  
        return
     end subroutine DumpGrid
                   
END Module FileIO