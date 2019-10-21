!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_sfcdiags

CONTAINS

  SUBROUTINE SFCDIAGS_RURAL(HFX, QFX, TSK, QSFC, CHS2, CQS2, T2, TH2, Q2,      &
      PSFC, CP, R_d, ROVCP,                            &
      ids, ide, jds, jde, kds, kde,                    &
      ims, ime, jms, jme, kms, kme,                    &
      its, ite, jts, jte, kts, kte,                    &
      XLAND) ! Add nakayoshi
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: HFX, QFX, TSK, QSFC, XLAND
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT) :: Q2, TH2, T2
    REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: PSFC, CHS2, CQS2
    REAL, INTENT(IN) :: CP, R_d, ROVCP

    ! LOCAL VARS
    INTEGER ::  i, j
    REAL    ::  rho

    do j = jts, jte
    do i = its, ite
    if (xland(i, j) < 1.5) then         ! add nakayoshi in open water
      rho = psfc(i, j) / (r_d * tsk(i, j)) ! density thermodynamics first law
      if (cqs2(i, j) < 1.e-5) then
        q2(i, j) = qsfc(i, j)
      else
        q2(i, j) = qsfc(i, j) - qfx(i, j) / (rho * cqs2(i, j))
      end if
      if (chs2(i, j) < 1.e-5) then
        t2(i, j) = tsk(i, j)
      else
        t2(i, j) = tsk(i, j) - hfx(i, j) / (rho * cp * chs2(i, j))
      end if
      th2(i, j) = t2(i, j) * (1.e5 / psfc(i, j)) ** rovcp
    end if
    end do
    end do

  END SUBROUTINE SFCDIAGS_RURAL

   SUBROUTINE SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2,TH2,Q2,       &
                     PSFC,CP,R_d,ROVCP,                            &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     )
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER,  INTENT(IN )   ::        ids,ide, jds,jde, kds,kde, &
                                        ims,ime, jms,jme, kms,kme, &
                                        its,ite, jts,jte, kts,kte
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN)                  ::                HFX, &
                                                              QFX, &
                                                              TSK, &
                                                             QSFC
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(INOUT)               ::                Q2, &
                                                             TH2, &
                                                              T2
      REAL,     DIMENSION( ims:ime, jms:jme )                    , &
                INTENT(IN)                  ::               PSFC, &
                                                             CHS2, &
                                                             CQS2
      REAL,     INTENT(IN   )               ::       CP,R_d,ROVCP
! LOCAL VARS
      INTEGER ::  I,J
      REAL    ::  RHO

      DO J=jts,jte
        DO I=its,ite
          RHO = PSFC(I,J)/(R_d * TSK(I,J)) !Density Thermodynamics first law
          if(CQS2(I,J).lt.1.E-5) then
             Q2(I,J)=QSFC(I,J)
          else
             Q2(I,J) = QSFC(I,J) - QFX(I,J)/(RHO*CQS2(I,J))
          endif
          if(CHS2(I,J).lt.1.E-5) then
             T2(I,J) = TSK(I,J)
          else
             T2(I,J) = TSK(I,J) - HFX(I,J)/(RHO*CP*CHS2(I,J))
          endif
          TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**ROVCP
        ENDDO
      ENDDO

  END SUBROUTINE SFCDIAGS

END MODULE module_sf_sfcdiags
