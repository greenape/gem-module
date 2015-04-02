
! Stand alone emulator
! using setup files produced by GEM-SA.exe
!
! Marc Kennedy
! m.kennedy@sheffield.ac.uk
! 26 Febuary 2004
! Modified to work as a python module by Jono Gray
! j.gray@soton.ac.uk
! 20/3/15

MODULE EMULATOR
real, allocatable :: xmod(:,:), ainv(:,:), g(:), ainvh(:,:), inv_hainvh(:,:)
real, allocatable :: h(:), betahat(:), roughness(:), min(:), max(:), t(:), ainvt(:)
real, allocatable :: tainvh(:), hdiff(:)
real prec, scalemean, scalesd
integer nmod, nreg, ninputs
END MODULE EMULATOR

SUBROUTINE INIT(path)
! open files, read data
USE EMULATOR

implicit none
character(*) path
integer i, j

IF(ALLOCATED(h)) DEALLOCATE(h, t, min, max, ainv, xmod, ainvh, betahat, roughness,&
g, inv_hainvh, ainvt, tainvh, hdiff)

!open emulator setup files
OPEN(UNIT=20, FILE=path//'emulator_ainv.txt', status='old')  ! inverse correlation matrix
OPEN(UNIT=21, FILE=path//'emulator_training_inputs.txt', status='old') ! training data inputs
OPEN(UNIT=24, FILE=path//'emulator_mu_out.dat', status='old') ! estimated regression parameters
OPEN(UNIT=25, FILE=path//'emulator_precision_out.dat', status='old') ! estimated GP prec = 1/variance
OPEN(UNIT=26, FILE=path//'emulator_ainvh.txt', status='old') ! ainv.H
OPEN(UNIT=27, FILE=path//'emulator_rough_out.dat', status='old') ! estimated function roughnesses
OPEN(UNIT=28, FILE=path//'emulator_g.txt', status='old') ! ainv.(y-Hb)
OPEN(UNIT=29, FILE=path//'emulator_inv_hainvh.txt', status='old') ! (H'.ainv.H)^{-1}
OPEN(UNIT=30, FILE=path//'emulator_minmax.txt', status='old') ! maxs and mins of each input
OPEN(UNIT=31, FILE=path//'emulator_scale.txt', status='old') ! mean and sd of output scaling

read(21,*)nmod, ninputs, nreg
ALLOCATE(h(nreg), t(nmod), min(ninputs), max(ninputs), ainv(nmod,nmod), xmod(nmod,ninputs), &
     ainvh(nmod,nreg), betahat(nreg), roughness(ninputs), g(nmod), inv_hainvh(nreg,nreg), &
     ainvt(nmod),tainvh(nreg), hdiff(nreg))

read(30,*)(min(i), i=1, ninputs)
read(30,*)(max(i), i=1, ninputs)
read(31,*)scalemean, scalesd

DO i = 1, nmod
  read(20,*)(ainv(i,j), j = 1, nmod)
  read(21,*)(xmod(i,j), j = 1, ninputs)
  read(26,*)(ainvh(i,j), j = 1, nreg)
ENDDO

read(24,*)(betahat(i), i=1, nreg)
read(25,*)prec
read(27,*)(roughness(i), i=1, ninputs)
read(28,*)(g(i), i=1, nmod)

DO i = 1, nreg
  read(29,*)(inv_hainvh(i, j), j = 1, nreg)
ENDDO

CLOSE(20)
CLOSE(21)
CLOSE(24)
CLOSE(25)
CLOSE(26)
CLOSE(27)
CLOSE(28)
CLOSE(29)
CLOSE(30)
CLOSE(31)

RETURN

END
!-------------------------------------------------------------------


SUBROUTINE NEXTPOINT(xpred, x, sd)
! reads next emulator prediction input and evaluates the emulator at that point
! mean and standard deviation returned through x and sd
USE EMULATOR
implicit none
integer i
real x, sd, aux, xpred(:)
                                                               

! read next point into xpred
! read(22,*)(xpred(i), i=1, ninputs)

! rescale inputs exactly as in the emulator building code
xpred = (xpred - min)/(max-min)

h(1) = 1.
IF (nreg.gt.1) THEN
  DO i=1, ninputs
    h(i+1) = xpred(i)  ! assumes linear regression terms only
  ENDDO
ENDIF

! compute t vector of correlations between xpred and training data inputs
DO i=1, nmod
  t(i) = exp(-SUM(roughness*(xpred-xmod(i,:))*(xpred-xmod(i,:))))
ENDDO

! compute mean
x = scalemean + scalesd*(DOT_PRODUCT(h, betahat) + DOT_PRODUCT(t, g))

! compute standard deviation
ainvt = MATMUL(ainv, t)
tainvh = MATMUL(t, ainvh)
hdiff = h - tainvh
aux = (1.0 - DOT_PRODUCT(t, ainvt) &
                  + DOT_PRODUCT(hdiff, MATMUL(inv_hainvh, hdiff)))/prec
IF (aux .le. 0.) THEN
  sd = 0.
ELSE
  sd = scalesd*SQRT(aux)
ENDIF
RETURN

END
!-------------------------------------------------------------------
