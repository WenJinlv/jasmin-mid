










c***********************************************************************
c
c    解Euler方程.
c
c***********************************************************************
      subroutine solvepde(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,dt,gamma,
     &  xo, yo, uo, dux, duy, flux, fluy, un)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
      double precision
     &    dt,gamma,DERMZY,PI
c
c variables in 2d indexed
c
      double precision
     &     xo(ifirst0-gcw0:ilast0+1+gcw0,
     &          ifirst1-gcw1:ilast1+1+gcw1),
     &     yo(ifirst0-gcw0:ilast0+1+gcw0,
     &          ifirst1-gcw1:ilast1+1+gcw1),
     &     uo(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,2+2),
     &     dux(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,2+2),
     &     duy(ifirst0-gcw0:ilast0+gcw0,
     &          ifirst1-gcw1:ilast1+gcw1,2+2),
     &     flux(ifirst0:ilast0+1,
     &          ifirst1:ilast1,2+2),
     &     fluy(ifirst0:ilast0,
     &          ifirst1:ilast1+1,2+2),
     &     un(ifirst0:ilast0,
     &          ifirst1:ilast1,2+2)
c
      integer J,K,L
      double precision
     &  DTX,DTY,ALENG,sztak,cztak,UU(4),
     &  RHOL,UXL,UYL,UNL,UTL,PREL,ALAML,
     &  XI2L,TEML0,TEML1,TEML2,TEML3,TGYLO,TGYL1,TGYL2,
     &  RHOR,UXR,UYR,UNR,UTR,PRER,ALAMR,
     &  XI2R,TEMR0,TEMR1,TEMR2,TEMR3,TGYRO,TGYR1,TGYR2,
     &  AREA,RR,CK
c
c***********************************************************************     
c
      PI = asin(1.0d0)*2.0
      CK = -2.0d0 + 2.0d0/(gamma-1.0d0)
C
C   X-DIRECTION FLUX
C
      DO K=ifirst1,ilast1
      DO J=ifirst0,ilast0+1
            DTX=XO(j,k+1)-XO(j,k  )
            DTY=YO(j,k+1)-YO(j,k  )
            ALENG=SQRT(DTX**2+DTY**2)
            sztak=DTY/ALENG
            cztak=DTX/ALENG     
          
            DO L=1,4
              UU(L)=UO(J-1,K,L)+0.5*(XO(J,K)-XO(J-1,K))*DUX(J-1,K,L)
            ENDDO
            RHOL=UU(1)
            UXL=UU(2)/UU(1)
            UYL=UU(3)/UU(1)
            UNL=sztak*UXL-cztak*UYL
            UTL=cztak*UXL+sztak*UYL           
            PREL=(gamma-1.)*(UU(4)-0.5*RHOL*(UXL*UXL+UYL*UYL))
            ALAML=0.5*RHOL/PREL
            XI2L=0.5*CK/ALAML
          
            TEML0=0.5*DERMZY(-SQRT(ALAML)*UNL)
            TEML1=UNL*TEML0+0.5*EXP(-ALAML*UNL*UNL)/SQRT(ALAML*PI)
            TEML2=UNL*TEML1+0.5*TEML0/ALAML
            TEML3=UNL*TEML2+1.0*TEML1/ALAML                            
           
            TGYLO=1.0
            TGYL1=UTL
            TGYL2=UTL*UTL+0.5/ALAML
          
            FLUX(J,K,1)=RHOL*TEML1  
            FLUX(J,K,2)=RHOL*( SZTAK*TEML2+CZTAK*TEML1*TGYL1)
            FLUX(J,K,3)=RHOL*(-CZTAK*TEML2+SZTAK*TEML1*TGYL1)          
            FLUX(J,K,4)=RHOL*(TEML3+TEML1*(TGYL2+XI2L))/2.
          
            DO L=1,4
               UU(L)=UO(J,K,L)-0.5*(XO(J,K)-XO(J-1,K))*DUX(J,K,L)
            ENDDO          
            RHOR=UU(1)
            UXR=UU(2)/UU(1)
            UYR=UU(3)/UU(1)
            UNR=sztak*UXR-cztak*UYR
            UTR=cztak*UXR+sztak*UYR               
            PRER=(gamma-1.)*(UU(4)-0.5*RHOR*(UXR*UXR+UYR*UYR))
            ALAMR=0.5*RHOR/PRER

            XI2R=0.5*CK/ALAMR
            TEMR0=0.5*DERMZY( SQRT(ALAMR)*UNR)
            TEMR1=UNR*TEMR0-0.5*EXP(-ALAMR*UNR*UNR)/SQRT(ALAMR*PI)
            TEMR2=UNR*TEMR1+0.5*TEMR0/ALAMR
            TEMR3=UNR*TEMR2+1.0*TEMR1/ALAMR

            TGYRO=1.0
            TGYR1=UTR
            TGYR2=UTR*UTR+0.5/ALAMR
                        
            FLUX(J,K,1)=FLUX(J,K,1)+RHOR*TEMR1
            FLUX(J,K,2)=FLUX(J,K,2)
     &                  +RHOR*( SZTAK*TEMR2+CZTAK*TEMR1*TGYR1) 
            FLUX(J,K,3)=FLUX(J,K,3)
     &                  +RHOR*(-CZTAK*TEMR2+SZTAK*TEMR1*TGYR1) 
            FLUX(J,K,4)=FLUX(J,K,4)
     &                  +RHOR*(TEMR3+TEMR1*(XI2R+TGYR2))/2
          
            DO L=1,4
               FLUX(J,K,L)=FLUX(J,K,L)*ALENG
            ENDDO

      ENDDO          
      ENDDO          
C
C  Y-DIRECTION FLUX
C
      DO K=ifirst1,ilast1+1
      DO J=ifirst0,ilast0

             DTX=XO(J  ,K)-XO(J+1,K)
             DTY=YO(J  ,K)-YO(J+1,K)
             ALENG=SQRT(DTX**2+DTY**2)

             sztak=DTY/ALENG
             cztak=DTX/ALENG
                
             DO L=1,4
                UU(L)=UO(J,K-1,L)+0.5*(YO(J,K)-YO(J,K-1))*DUY(J,K-1,L)
             ENDDO

             RHOL=UU(1)
             UXL=UU(2)/UU(1)
             UYL=UU(3)/UU(1)
             UNL=sztak*UXL-cztak*UYL
             UTL=cztak*UXL+sztak*UYL           
             PREL=(gamma-1.)*(UU(4)-0.5*RHOL*(UXL*UXL+UYL*UYL))
             ALAML=0.5*RHOL/PREL
             XI2L=0.5*CK/ALAML
          
             TEML0=0.5*DERMZY(-SQRT(ALAML)*UNL)
             TEML1=UNL*TEML0+0.5*EXP(-ALAML*UNL*UNL)/SQRT(ALAML*PI)
             TEML2=UNL*TEML1+0.5*TEML0/ALAML
             TEML3=UNL*TEML2+1.0*TEML1/ALAML                            
          
             TGYLO=1.0
             TGYL1=UTL
             TGYL2=UTL*UTL+0.5/ALAML
          
             FLUY(J,K,1)=RHOL*TEML1  
             FLUY(J,K,2)=RHOL*( SZTAK*TEML2+CZTAK*TEML1*TGYL1)
             FLUY(J,K,3)=RHOL*(-CZTAK*TEML2+SZTAK*TEML1*TGYL1)          
             FLUY(J,K,4)=RHOL*(TEML3+TEML1*(TGYL2+XI2L))/2.          

             DO L=1,4
                UU(L)=UO(J,K,L)-0.5*(YO(J,K)-YO(J,K-1))*DUY(J,K,L)
             ENDDO          
             RHOR=UU(1)
             UXR=UU(2)/UU(1)
             UYR=UU(3)/UU(1)
             UNR=sztak*UXR-cztak*UYR
             UTR=cztak*UXR+sztak*UYR               
             PRER=(gamma-1.)*(UU(4)-0.5*RHOR*(UXR*UXR+UYR*UYR))
             ALAMR=0.5*RHOR/PRER

             XI2R=0.5*CK/ALAMR
             TEMR0=0.5*DERMZY( SQRT(ALAMR)*UNR)
             TEMR1=UNR*TEMR0-0.5*EXP(-ALAMR*UNR*UNR)/SQRT(ALAMR*PI)
             TEMR2=UNR*TEMR1+0.5*TEMR0/ALAMR
             TEMR3=UNR*TEMR2+1.0*TEMR1/ALAMR

             TGYRO=1.0
             TGYR1=UTR
             TGYR2=UTR*UTR+0.5/ALAMR

             FLUY(J,K,1)=FLUY(J,K,1)+RHOR*TEMR1
             FLUY(J,K,2)=FLUY(J,K,2)
     &                   +RHOR*( SZTAK*TEMR2+CZTAK*TEMR1*TGYR1) 
             FLUY(J,K,3)=FLUY(J,K,3)
     &                   +RHOR*(-CZTAK*TEMR2+SZTAK*TEMR1*TGYR1) 
             FLUY(J,K,4)=FLUY(J,K,4)
     &                   +RHOR*(TEMR3+TEMR1*(XI2R+TGYR2))/2

             DO L=1,4
                FLUY(J,K,L)=FLUY(J,K,L)*ALENG
             ENDDO

      ENDDO 
      ENDDO 
c
      DO K=ifirst1,ilast1
      DO J=ifirst0,ilast0
          AREA=0.5*((XO(j,k)-XO(j+1,k+1))*(YO(j+1,k)-YO(j,k+1))
     &             -(YO(j,k)-YO(j+1,k+1))*(XO(j+1,k)-XO(j,k+1)))          
          RR=DT/AREA
          DO L=1,4
              UN(J,K,L)=UO(J,K,L)-RR*(FLUX(J+1,K,L)-FLUX(J,K,L)) 
     &                           -RR*(FLUY(J,K+1,L)-FLUY(J,K,L))             
          ENDDO
      ENDDO
      ENDDO
c
      return
      end
c
      FUNCTION DERMZY(X)
      implicit real*8 (a-h,o-z)
      Z=ABS(X)
      T=1.0d0/(1.0d0+0.5d0*Z)
      DERMZY = T*EXP(-Z*Z-1.26551223d0+
     &        T*(1.000023680d0+T*(0.374091960d0+
     &        T*(0.096784180d0+T*(-0.186288060d0
     &       +T*(.278868070d0+T*(-1.135203980d0+
     &        T*(1.48851870d0+T*(-.822152230d0
     &       +T*.170872770d0)))))))))
      IF(X.LT.0.0) DERMZY=2.0d0-DERMZY
c
      RETURN
      END


