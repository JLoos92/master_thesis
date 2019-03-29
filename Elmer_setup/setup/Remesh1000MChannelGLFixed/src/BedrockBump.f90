FUNCTION BedrockBump( Model, nodenumber, Time) RESULT(BedBump)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: BedSol
   INTEGER, POINTER :: BedPerm(:)
   REAL(kind=dp), POINTER :: BedVal(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: BedBump, Time, x,y,y0,x0,sigmax,sigmay,Ampl,Bump
   REAL(KIND=dp) :: MaxGrowTime
   LOGICAL :: FirstTime=.True., UnFoundFatal

   SAVE FirstTime
   MaxGrowTime = 50.0 !time it takes for the bump to grow
   Ampl = 150.0
   sigmax = 500
   sigmay = 500
   x0 = 1056000
   y0 = 0
   x = Model % Nodes % x(nodenumber) ! get coordinates
   y = Model % Nodes % y(nodenumber) ! get coordinates
   
   !Get Bedrock without catching any error messages if fields don't exist
   BedSol => VariableGet( Model % Variables, 'Bedrock',UnFoundFatal=UnFoundFatal)
   BedPerm => BedSol % Perm
   BedVal => BedSol % Values

   !Bump = exp(-((x-x0)**2/(2*sigmax**2)))
   !NORMAL BUMP
   if (Time <= MaxGrowTime) then
      Bump =(Time/MaxGrowTime)*Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   else
      Bump =Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   end if
   ! Grow and decay bump
   !if (Time <= MaxGrowTime.AND.Time < DecayStart) then
      !Bump =(Time/MaxGrowTime)*Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   !elseif (Time > MaxGrowTime.AND.Time < DecayStart) then
      !Bump =Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   !elseif (Time >= DecayStart.AND. Time < DecayStart+MaxGrowTime) then
      !Bump =((Time-DecayStart+MaxGrowTime)/(-MaxGrowTime))*Ampl * exp(-((x-x0)**2/(2*sigmax**2) + (y-y0)**2/(2*sigmay**2))) 
   !else 
      !Bump=0
   !end if
   !write(*,*) 'BUMP',Bump
   BedBump = BedVal(BedPerm(nodenumber))+Bump

END FUNCTION BedrockBump
