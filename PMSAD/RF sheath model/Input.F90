Module Input
	use constants
	implicit none
           Integer(4),parameter :: Dims=2
           Integer(4),parameter :: InputNx=61,InputNy=151
           Real(8),parameter,private :: XLength=0.005d0
           Real(8),parameter :: Inputdx = XLength/dble(InputNx-1)
           Real(8),parameter :: Inputdt=1.d-11
           Integer(4),parameter :: InputGap=0
           Integer(4),parameter :: Nelech=1   
           Integer(4),parameter :: ParticlePerGrid=200   
           Real(8),parameter :: InitDensity=1.d18        
           Real(8) :: BxField(InputNx,InputNy)=0.d0
           Real(8) :: ByField(InputNx,InputNy)=0.d0
           Real(8) :: BzField(InputNx,InputNy)=0.d0
           Integer(4) :: RFKind=21
           Real(8) :: RFFrequency(2)=(/37.d6,60.d6/),RfVol(2)=(/300.d0,0.d0/)
           Real(8) :: YFrequency=60.d6,YVol=0.d0 
           Integer(4) :: InputNGas=1
           Integer(4) :: GasDefinition(1:3)=(/1,0,0/)
           Real(8) :: GasPressure(1:3)=(/30.d0,0.d0,0.d0/),GasTemperature(1:3)=(/1160500.d0,300.d0,300.d0/)			
end Module Input




 
