Module ParticleModule
    Use TypeModule
    Use Constants
    Implicit none
    Contains
	
    Subroutine  AddParticle(OneParticle,InputParticle)
        Implicit none
        Type(ParticleOne),intent(in) ::  OneParticle
        Type(ParticleBundle),intent(out) ::  InputParticle
        InputParticle%NPar=InputParticle%NPar+1
        InputParticle%PhaseSpace(InputParticle%NPar)=OneParticle
        return   
    End  Subroutine  AddParticle
     
    Subroutine  DelParticle(NDel,InputParticle)
        Implicit none
        integer(4),intent(in) ::  NDel
        Type(ParticleBundle),intent(out) ::  InputParticle
        InputParticle%PhaseSpace(NDel)=InputParticle%PhaseSpace(InputParticle%NPar)
        InputParticle%NPar=InputParticle%NPar-1
        return   
    End  Subroutine  DelParticle
             
    Subroutine  UpdatePositionBundle(InputParticle,Sign)
        Implicit none
        Type(ParticleBundle), intent(out) ::  InputParticle
        Integer(4),intent(in) :: Sign
        Real(8) :: XFactor
        Integer(4) :: i
        Select case (Sign)
            Case(1_4)
                XFactor=InputParticle%XFactor
            Case(-1_4)
                 XFactor=1.d0/InputParticle%XFactor
        End Select       
        do i=1,InputParticle%NPar
            InputParticle%PhaseSpace(i)%X=XFactor*InputParticle%PhaseSpace(i)%X
            InputParticle%PhaseSpace(i)%Y=XFactor*InputParticle%PhaseSpace(i)%Y
        end do 
        return 
    End Subroutine  UpdatePositionBundle
             
    Subroutine  UpdateVelocityBundle(InputParticle,Sign)
        Implicit none
        Type(ParticleBundle), intent(out) ::  InputParticle
        Integer(4),intent(in) :: Sign
        Real(8) :: VFactor
        Integer(4) :: i
        Select case (Sign)
            Case(1_4)
                VFactor=InputParticle%VFactor
            Case(-1_4)
                VFactor=1.d0/InputParticle%VFactor
        End Select 
        do i=1,InputParticle%NPar
            InputParticle%PhaseSpace(i)%Vx=VFactor*InputParticle%PhaseSpace(i)%Vx
            InputParticle%PhaseSpace(i)%Vy=VFactor*InputParticle%PhaseSpace(i)%Vy
            InputParticle%PhaseSpace(i)%Vz=VFactor*InputParticle%PhaseSpace(i)%Vz
        end do 
        return 
	End Subroutine  UpdateVelocityBundle      
                                                                     
           Subroutine CalEnergy(Mass,InputParticle,Energy) 
              Implicit none
              Real(8),intent(in):: Mass 
              Type(ParticleOne),intent(in) :: InputParticle
              Real(8),intent(out) ::Energy
              Energy=0.5d0*Mass*(InputParticle%Vx*InputParticle%Vx+InputParticle%Vy*InputParticle%Vy +InputParticle%Vz*InputParticle%Vz)/JtoeV
             return 
           end  Subroutine  CalEnergy
          
           
End  Module ParticleModule

