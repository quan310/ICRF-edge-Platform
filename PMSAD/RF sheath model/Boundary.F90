Module ParticleBoundary
    Use ParticleModule
    implicit none
    contains
	
    Subroutine  AborptionXYInit(InputXmin,InputXmax,InputYmin,InputYmax,Gap)
                Implicit none
                Real(8),intent(in) :: InputXmin,InputXmax,InputYmin,InputYmax,charge,dx
                
                Type(ParticleBundle),intent(inout) :: InputParticle
                Real(8),save :: Xmin,Xmax,Ymin,Ymax,Xbottomgap,Xleftgap
                Integer(4):: i,Gap
                Real(8),intent(inout) ::  QBottom,QTop
                Xmin=InputXmin
                Xmax=InputXmax
                Ymin=InputYmin
                Ymax=InputYmax
                Xbottomgap=Ymax-dble(Gap)
                Xleftgap=dble(Gap)
                return
                Entry AborptionXY(InputParticle,charge,QBottom,QTop,dx)
                 QBottom=0.d0
                 QTop=0.d0
                 do i=InputParticle%NPar,1,-1
        
                       If(InputParticle%PhaseSpace(i)%Y>Ymax) then
                              Call DelParticle(i,InputParticle)
                             
                       else if (InputParticle%PhaseSpace(i)%Y<Ymin) then
                           
				InputParticle%PhaseSpace(i)%Y=-InputParticle%PhaseSpace(i)%Y+1.d-12
				InputParticle%PhaseSpace(i)%Vy=-(InputParticle%PhaseSpace(i)%VY)-1.d-12
                                  If  (InputParticle%PhaseSpace(i)%X<Xmin) then
                                             if(InputParticle%PhaseSpace(i)%y<29.0) then
                                           	  InputParticle%PhaseSpace(i)%X=-InputParticle%PhaseSpace(i)%X+1.d-12
				                  InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                                             else  if(InputParticle%PhaseSpace(i)%y<50.0) then
                                                      Call DelParticle(i,InputParticle)  
                                             else 
                                                QBottom=QBottom+1.d0            
                                                Call DelParticle(i,InputParticle)  
                                             end if

                                  else If  (InputParticle%PhaseSpace(i)%X>Xmax) then

    			          	InputParticle%PhaseSpace(i)%X=2*Xmax-InputParticle%PhaseSpace(i)%X-1.d-12
			        	InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                     
                                 End if
     


                       else If  (InputParticle%PhaseSpace(i)%X<Xmin) then
                                  if(InputParticle%PhaseSpace(i)%y<29.0) then
                                         	InputParticle%PhaseSpace(i)%X=-InputParticle%PhaseSpace(i)%X+1.d-12
				               InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                                   else  if(InputParticle%PhaseSpace(i)%y<50.0) then
                                             Call DelParticle(i,InputParticle)  
                                 else 
                                                QBottom=QBottom+1.d0            
                                                Call DelParticle(i,InputParticle)  
                                  end if
                              If  (InputParticle%PhaseSpace(i)%Y>50.0) then

                                 QBottom=QBottom+1.d0 
                              end if          
                              Call DelParticle(i,InputParticle)  
              
                       else If  (InputParticle%PhaseSpace(i)%X>Xmax) then

   
				InputParticle%PhaseSpace(i)%X=2*Xmax-InputParticle%PhaseSpace(i)%X-1.d-12
				InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                       End if
                 End do 


                    return
          End  Subroutine AborptionXYInit

	Subroutine  AborptionXYInit1(InputXmin,InputXmax,InputYmin,InputYmax,Gap)        
		use input
        Implicit none		
        Real(8),intent(in) :: InputXmin,InputXmax,InputYmin,InputYmax,charge,dx
        Type(ParticleBundle),intent(inout) :: InputParticle
        Real(8),save :: Xmin,Xmax,Ymin,Ymax,Xbottomgap,Xleftgap
        Integer(4) :: rhoagain(inputNy)
        Integer(4):: i,Gap
        Real(8),intent(inout) ::  QBottom,QTop		
        Xmin=InputXmin
        Xmax=InputXmax
        Ymin=InputYmin
        Ymax=InputYmax
        Xbottomgap=Ymax-dble(Gap)
        Xleftgap=dble(Gap)		
        return	
		
		Entry AborptionXY1(InputParticle,charge,QBottom,QTop,dx,rhoagain)          
			QBottom=0.d0
            QTop=0.d0
            rhoagain=0.0           
			do i = InputParticle%NPar,1,-1     
                If(InputParticle%PhaseSpace(i)%Y>Ymax) then
                    Call DelParticle(i,InputParticle)
                else if (InputParticle%PhaseSpace(i)%Y<Ymin) then                        
				InputParticle%PhaseSpace(i)%Y=-InputParticle%PhaseSpace(i)%Y !Ymax+InputParticle%PhaseSpace(i)%Y !
				InputParticle%PhaseSpace(i)%Vy=-(InputParticle%PhaseSpace(i)%VY)-1.d-12                
				end if     			
			End do 

            do i=InputParticle%NPar,1,-1       
                If(InputParticle%PhaseSpace(i)%Y>Ymax) then
                    Call DelParticle(i,InputParticle)
                else if (InputParticle%PhaseSpace(i)%Y<Ymin) then                          
					InputParticle%PhaseSpace(i)%Y=-InputParticle%PhaseSpace(i)%Y!Ymax+InputParticle%PhaseSpace(i)%Y !
					inputParticle%PhaseSpace(i)%Vy=-(InputParticle%PhaseSpace(i)%VY)-1.d-12
                end if     
                If  (InputParticle%PhaseSpace(i)%X<Xmin) then
                    if(InputParticle%PhaseSpace(i)%y<30.0) then                       
					InputParticle%PhaseSpace(i)%X=-InputParticle%PhaseSpace(i)%X+1.d-12
				    InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12                    
					else                        
					QBottom=QBottom+1.d0   
					Call DelParticle(i,InputParticle) 
					rhoagain(ceiling(InputParticle%PhaseSpace(i)%Y))=rhoagain(ceiling(InputParticle%PhaseSpace(i)%Y))+1
					end if
					else If  (InputParticle%PhaseSpace(i)%X>Xmax) then
					InputParticle%PhaseSpace(i)%X = 2*Xmax - InputParticle%PhaseSpace(i)%X-1.d-12
					InputParticle%PhaseSpace(i)%Vx = -(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                End if
            End do 

			do i=InputParticle%NPar,1,-1      
                if(InputParticle%PhaseSpace(i)%Y>Ymax.or.InputParticle%PhaseSpace(i)%Y<Ymin.or.InputParticle%PhaseSpace(i)%X<Xmin.or.InputParticle%PhaseSpace(i)%X>Xmax) then
					Call DelParticle(i,InputParticle)  
				end if 
			End do 
		return
    End  Subroutine AborptionXYInit1

	Subroutine  AborptionXYInit2(InputXmin,InputXmax,InputYmin,InputYmax,Gap)
        Implicit none
        Real(8),intent(in) :: InputXmin,InputXmax,InputYmin,InputYmax,charge,dx              
        Type(ParticleBundle),intent(inout) :: InputParticle
        Real(8),save :: Xmin,Xmax,Ymin,Ymax,Xbottomgap,Xleftgap
        Integer(4):: i,Gap
        Real(8),intent(inout) ::  QBottom,QTop
                Xmin=InputXmin
                Xmax=InputXmax
                Ymin=InputYmin
                Ymax=InputYmax
                Xbottomgap=Ymax-dble(Gap)
                Xleftgap=dble(Gap)
        return
		Entry AborptionXY2(InputParticle,charge,QBottom,QTop,dx)
            QBottom=0.d0
            QTop=0.d0
            do i=InputParticle%NPar,1,-1        
                     If(InputParticle%PhaseSpace(i)%Y>Ymax) then
                                 Call DelParticle(i,InputParticle)
                       else if (InputParticle%PhaseSpace(i)%Y<Ymin) then                           
				InputParticle%PhaseSpace(i)%Y=-InputParticle%PhaseSpace(i)%Y				!Ymax+InputParticle%PhaseSpace(i)%Y !
				InputParticle%PhaseSpace(i)%Vy=-(InputParticle%PhaseSpace(i)%VY)-1.d-12
                        end if     
                  End do 
               do i=InputParticle%NPar,1,-1       
                        If(InputParticle%PhaseSpace(i)%Y>Ymax) then
                                Call DelParticle(i,InputParticle)
                       else if (InputParticle%PhaseSpace(i)%Y<Ymin) then                          
				InputParticle%PhaseSpace(i)%Y=-InputParticle%PhaseSpace(i)%Y!Ymax+InputParticle%PhaseSpace(i)%Y !
				inputParticle%PhaseSpace(i)%Vy=-(InputParticle%PhaseSpace(i)%VY)-1.d-12
                        end if     
                      If  (InputParticle%PhaseSpace(i)%X<Xmin) then
                                 if(InputParticle%PhaseSpace(i)%y<30.0) then
                                         	InputParticle%PhaseSpace(i)%X=-InputParticle%PhaseSpace(i)%X+1.d-12
				               InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                           else 
                                 QBottom=QBottom+1.d0   
                                 Call DelParticle(i,InputParticle)  
                             end if
                       else If  (InputParticle%PhaseSpace(i)%X>Xmax) then
				InputParticle%PhaseSpace(i)%X=2*Xmax-InputParticle%PhaseSpace(i)%X-1.d-12
				InputParticle%PhaseSpace(i)%Vx=-(InputParticle%PhaseSpace(i)%Vx)-1.d-12
                       End if
                 End do 
				do i=InputParticle%NPar,1,-1        
                   if(InputParticle%PhaseSpace(i)%Y>Ymax.or.InputParticle%PhaseSpace(i)%Y<Ymin.or.InputParticle%PhaseSpace(i)%X<Xmin.or.InputParticle%PhaseSpace(i)%X>Xmax) then
					Call DelParticle(i,InputParticle)  
					end if 
				end do
        return
   End  Subroutine AborptionXYInit2
   
End Module ParticleBoundary

Module RFSource
    Use Constants 
    implicit none
    contains
    Subroutine RfVoltageXInitilalization(InputStyle,InputFrequency,InputVoltage,Inputdt,Period)
               Implicit none
               Integer(4), intent(in) :: InputStyle
               Real(8),intent(in) :: InputFrequency(2),InputVoltage(2),Inputdt
               Integer(4), intent(out) ::  Period
               Real(8),intent(out) :: Top,Bottom  
                
               Integer(4),save :: Style,Timer=0
               Real(8),save :: Frequency(2),Voltage(2),dt,FMin
               Style=InputStyle
               Frequency=InputFrequency
               Voltage=InputVoltage
               dt=Inputdt
               !FMin=MINVAL(Frequency)
               Period=Int(1.d0/Frequency(1)/dt)
               return
                
       Entry  RfVoltageX(Bottom,Top)
               Select case (Style)
                      case (11)
                               Top=0.d0 
                               Bottom=Voltage(1)*DSin(2*PI*Frequency(1)*dt*Dble(Timer))
                      case (21)
                               Top=0.d0 
                               Bottom=Voltage(1)*DSin(2*PI*Frequency(1)*dt*Dble(Timer))+Voltage(2)*DSin(2*PI*Frequency(2)*dt*Dble(Timer))
                      case(22)
                               Top=Voltage(2)*DSin(2*PI*Frequency(2)*dt*Dble(Timer))
                               Bottom=Voltage(1)*DSin(2*PI*Frequency(1)*dt*Dble(Timer))  
              end  Select
               Timer=Timer+1
               return
       Entry  RfVoltageX2(Bottom,Top)
               Select case (Style)
                      case (11)
                               Top=0.d0 
                               Bottom=Voltage(1)*DSin(2*PI*Frequency(1)*dt*Dble(Timer))
                      case (21)
                               Top=0.d0 
                               Bottom=Voltage(1)*DSin(2*PI*Frequency(1)*dt*Dble(Timer))+Voltage(2)*DSin(2*PI*Frequency(2)*dt*Dble(Timer))
                      case(22)
                               Top=Voltage(2)*DSin(2*PI*Frequency(2)*dt*Dble(Timer))
                               Bottom=Voltage(1)*DSin(2*PI*Frequency(1)*dt*Dble(Timer))  
              end  Select
              
               return
    end subroutine  RfVoltageXInitilalization
   
    Subroutine RfVoltageYInitilalization(InputFrequency,InputVoltage,Inputdt)
               Implicit none
               Real(8),intent(in) :: InputFrequency,InputVoltage,Inputdt
               Real(8),intent(out) :: Left,Right  
               Integer(4),save :: Timer=0
               Real(8),save :: Frequency,Voltage,dt
               Frequency=InputFrequency
               Voltage=InputVoltage
               dt=Inputdt
               return
         Entry  RfVoltageY(Left,Right)
                     Left=0
                     Right=Voltage*DSin(2*PI*Frequency*dt*Dble(Timer))
                     Timer=Timer+1
                     
               return
    end subroutine  RfVoltageYInitilalization
end  Module RFSource

