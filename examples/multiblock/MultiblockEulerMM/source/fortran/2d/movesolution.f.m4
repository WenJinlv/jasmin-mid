define(NDIM,2)dnl
define(REAL,`double precision')dnl
include(JASMIN_FORTDIR/pdat_m4arrdim2d.i)dnl

c***********************************************************************
c
c    移动网格后, 更新数值解.
c
c***********************************************************************
      subroutine movesolution(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  gcw0,gcw1,
     &  xo, yo, xn,yn, uo, un, dux, duy, flux, fluy)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1
      integer gcw0,gcw1
c
c variables in 2d indexed
c
      REAL
     &     xo(NODE2dVECG(ifirst,ilast,gcw)),
     &     yo(NODE2dVECG(ifirst,ilast,gcw)),
     &     xn(NODE2dVECG(ifirst,ilast,gcw)),
     &     yn(NODE2dVECG(ifirst,ilast,gcw)),
     &     uo(CELL2dVECG(ifirst,ilast,gcw),NDIM+2),
     &     dux(CELL2dVECG(ifirst,ilast,gcw),NDIM+2),
     &     duy(CELL2dVECG(ifirst,ilast,gcw),NDIM+2),
     &     flux(SIDE2d0(ifirst,ilast,0),NDIM+2),
     &     fluy(SIDE2d1(ifirst,ilast,0),NDIM+2),
     &     un(CELL2d(ifirst,ilast,0),NDIM+2)
c
      integer J,K,L
      REAL
     &  DTX,DTY,ALENG,sztak,cztak,SPX,SPY,SPN,UUL,UUR,
     &  AREAO,AREAN,UNN
c
c***********************************************************************     
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
          SPX=.5*(XO(J,K)+XO(J,K+1)-XN(J,K)-XN(J,K+1))  
          SPY=.5*(YO(J,K)+YO(J,K+1)-YN(J,K)-YN(J,K+1))          
          SPN=sztak*SPX-cztak*SPY                     
          DO L=1,4
             UUL=UO(J-1,K,L)+0.5*(XO(J,K)-XO(J-1,K))*DUX(J-1,K,L)
             UUR=UO(J,K,L)-0.5*(XO(J+1,K)-XO(J,K))*DUX(J,K,L)          
             FLUX(J,K,L)=.5*(SPN*(UUR+UUL)-ABS(SPN)*(UUR-UUL))*ALENG
          ENDDO
      ENDDO          
      ENDDO          
C
C  Y-DIRECTION FLUX
C
      DO K=ifirst1,ilast1+1
      DO J=ifirst0,ilast0
          DTX=XO(j,k)-XO(j+1,k)
          DTY=YO(j,k)-YO(j+1,k)
          ALENG=SQRT(DTX**2+DTY**2)
          sztak=DTY/ALENG
          cztak=DTX/ALENG    
          SPX=.5*(XO(J,K)+XO(J+1,K)-XN(J,K)-XN(J+1,K))  
          SPY=.5*(YO(J,K)+YO(J+1,K)-YN(J,K)-YN(J+1,K))          
          SPN=sztak*SPX-cztak*SPY                     
          DO L=1,4
             UUL=UO(J,K-1,L)+0.5*(YO(J,K)-YO(J,K-1))*DUY(J,K-1,L)
             UUR=UO(J,K,L)-0.5*(YO(J,K+1)-YO(J,K))*DUY(J,K,L)          
             FLUY(J,K,L)=.5*(SPN*(UUR+UUL)-ABS(SPN)*(UUR-UUL))*ALENG
          ENDDO            
      ENDDO
      ENDDO
C
      DO K=ifirst1,ilast1
      DO J=ifirst0,ilast0
          AREAO=0.5*((XO(j,k)-XO(j+1,k+1))*(YO(j+1,k)-YO(j,k+1))
     &             -(YO(j,k)-YO(j+1,k+1))*(XO(j+1,k)-XO(j,k+1)))
          AREAN=0.5*((XN(j,k)-XN(j+1,k+1))*(YN(j+1,k)-YN(j,k+1))
     &             - (YN(j,k)-YN(j+1,k+1))*(XN(j+1,k)-XN(j,k+1)))               
          DO L=1,4        
             UNN=UO(J,K,L)*AREAO-(FLUX(J+1,K,L)-FLUX(J,K,L))
     &                          -(FLUY(J,K+1,L)-FLUY(J,K,L))
             UN(J,K,L)=UNN/AREAN                                
          ENDDO
      ENDDO
      ENDDO
c
      return
      end

