define(riemann_solve,`dnl
      do ic$3=ifirst$3-$6,ilast$3+$6
         do ic$2=ifirst$2-$5,ilast$2+$5
           do ie$1=ifirst$1,ilast$1+1

             advecspeed= (xlo($2)+dx($2)*(ic$2-ifirst$2+0.5))*omega($3)
     &                  -(xlo($3)+dx($3)*(ic$3-ifirst$3+0.5))*omega($2)
 
             if (advecspeed.ge.zero) then
               riemst= trlft$1(ie$1,$4)
             else
               riemst= trrgt$1(ie$1,$4)
             endif
             flux$1(ie$1,$4)= dt*riemst*advecspeed

           enddo
         enddo
      enddo
')dnl
define(correc_flux2d,`dnl
c   correct the $1-direction with $3-fluxes
      do ic$5=ifirst$5-1,ilast$5+1
         do ic$3=ifirst$3,ilast$3
           do ic$1=ifirst$1-1,ilast$1+1
             trnsvers=
     &           (flux$3(ic$3+1,$4)-flux$3(ic$3,$4))*half/dx($3)
c
             ttracelft$1(ic$1+1,$2)=tracelft$1(ic$1+1,$2)
     &                                    - trnsvers
             ttracergt$1(ic$1  ,$2)=tracergt$1(ic$1  ,$2)
     &                                    - trnsvers
           enddo
         enddo
      enddo
')dnl
define(correc_flux3d,`dnl
c   correct the $1-direction with $2$3-fluxes
      do ic$1=ifirst$1-1,ilast$1+1
         do ic$3=ifirst$3,ilast$3
           do ic$2=ifirst$2,ilast$2
             trnsvers=half*(
     &          (flux$4(ic$2+1,$6)-flux$4(ic$2,$6))/dx($2)+
     &          (flux$5(ic$3+1,$7)-flux$5(ic$3,$7))/dx($3))
  
             tracelft$1(ic$1+1,ic$2,ic$3)=tracelft$1(ic$1+1,ic$2,ic$3)
     &                                           - trnsvers
             tracergt$1(ic$1  ,ic$2,ic$3)=tracergt$1(ic$1  ,ic$2,ic$3)
     &                                           - trnsvers
           enddo
         enddo
      enddo
')dnl
