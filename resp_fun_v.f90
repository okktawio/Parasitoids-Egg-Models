real*8 function respuesta_funcional_1(Nt, a)

!Respuesta funcional 1
  
  implicit none
  real*8 a, Nt

  respuesta_funcional_1  = a * Nt
  return
end function respuesta_funcional_1

real*8 function respuesta_funcional_2(Nt, As, a)

!Respuesta funcional 2
  
  implicit none
  real*8 a, Nt, As
  
  respuesta_funcional_2  = (a * Nt)/(1 + a * Nt / As)
  return
end function respuesta_funcional_2


real*8 function respuesta_funcional_3(Nt, As, b, m)

!Respuesta funcional 3
  
  implicit none
  real*8 a, Nt, b, m, As
  
  a =m * Nt + b;
  if(a > 1.0) a = 1.0
  respuesta_funcional_3  = (a * Nt)/(1 + a * Nt / As)
  return
end function respuesta_funcional_3

real*8 function respuesta_funcional_4(Nt, As, b, m, cuh)

!Respuesta funcional 3 con aprendizaje acumulativo
  
  implicit none
  real*8 a, Nt, b, m, As, cuh
  
  a =m * cuh + b;
  if(a > 1.0) a = 1.0
  respuesta_funcional_4 = (a * Nt)/(1 + a * Nt / As)
  return
end function respuesta_funcional_4


real*8 function binomial(x, n, p)

!c binomial log-likelihood function     

!c  updated 17/01/2007. dh. 

  implicit none
  real*8 like, p
  integer x,n
!  logical not_scalar_n,not_scalar_p
  real*8 factln
  real*8 infinity
  parameter (infinity = 1.7976931348623157d300)
  
  like = 0.0

!  print *,"binomial 0",x,p,n,p/n

  p = p / n
  
!  print *,"binomial 1",x,p,n
  
  if ((x < 0) .or. (n < 0) .or. (x > n)) then
     like = -infinity
  else
     if (p < 0.00000001) then
        p = 0.00000001
     else if(p > 0.99999999) then
        p = 0.99999999
     end if
     like = x*log(p) + (n-x)*log(1.-p)
     like = like + factln(n) - factln(x) - factln(n-x) 
  end if
  binomial = like
  return
end function binomial


real*8 function poisson(x, mu)

! Poisson log-likelihood function      
! UPDATED 1/16/07 AP

  implicit none
  integer x
  real*8 like, mu, infinity, factln
  parameter (infinity = 1.7976931348623157d300)
  
  if (mu < 0.0) then
     like = -infinity
  else if (x < 0.0) then
     like = -infinity
  else if (.not.((x == 0.0) .and. (mu == 0.0))) then
     like = (x*log(mu) - mu) - factln(x)
  else 
     like = 0.0
  endif
  !poisson = exp(like)
  poisson = like
  return
end function poisson

real*8 function factln(n) 
  !c uses gammln returns ln(n!). 
  
  integer n 
  real*8 a(100),gammln, pass_val 
  real*8 infinity
  parameter (infinity = 1.7976931348623157d300)
  
  save a 
  !c initialize the table to negative values. 
  data a/100*-1./ 
  pass_val = n + 1
  if (n.lt.0) then
     !c        write (*,*) 'negative factorial in factln' 
     factln=-infinity
     return
  endif
  !c in range of the table. 
  if (n.le.99) then
     !c if not already in the table, put it in.
     if (a(n+1).lt.0.) a(n+1)=gammln(pass_val) 
     factln=a(n+1) 
  else 
     !c out of range of the table. 
     factln=gammln(pass_val) 
  endif
  return 
end function factln

real*8 function clip(x)
  real*8 x, tolmin, tolmax
  parameter (tolmin = 1.0d-10)
  parameter (tolmax = 1-(1.0d-10))
  
  if (x > tolmax) then 
     clip = tolmax
     return
  end if
  if (x < tolmin) then
     clip = tolmin
     return
  end if
  clip = x
  return
end function clip

real*8 function logit(x)
  real*8 x, p, clip
  p = clip(x)
  logit = log(p/(1-p))
  return
end function logit

real*8 function normal(x,mu,tau)
  ! Normal log-likelihood function      
  real*8 x, mu, tau, like
  !REAL*8 mu_tmp, tau_tmp
  real*8 pi
  parameter (pi=3.141592653589793238462643d0) 
  real*8 infinity
  parameter (infinity = 1.7976931348623157d300)
  !!print *,"x",x,"mu",mu,"tau",tau
  like = 0.0
  if (tau <= 0.0) then
     like = -infinity
     return
  end if
  like = like - 0.5 * tau * (x-mu)**2
  like = like + 0.5 * log(0.5 * tau/pi)

  !normal = exp(like)
  normal = like
  !!print *,exp(like)
  return
end function normal


real*8 function gammln(xx) 
  !c returns the value ln[gamma(xx)] for xx > 0. 
  
  real*8 xx
  integer j 
  real*8 ser, tmp,x,y
  real*8 stp, cof(6)
  
  !c internal arithmetic will be done in real*8, 
  !c a nicety that you can omit if five-figure accuracy is good enough. 
  save cof,stp

  cof = (/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,&
       -1.231739572450155d0,0.1208650973866179d-2,-.5395239384953d-5/)
  stp = 2.5066282746310005d0
  x=xx
  y=x 
  tmp=x+5.5d0 
  tmp=(x+0.5d0)*log(tmp)-tmp 
  ser=1.000000000190015d0 
  do j=1,6
     y=y+1.d0 
     ser=ser+cof(j)/y 
  enddo
  gammln=tmp+log(stp*ser/x) 
  return 
end function gammln

module modelos_bichos
!modulo con los modelos de bichos
!1 limitado por huevos sin produccion
!2 limitado por huevos con produccion constante
!3 limitado por huevos con produccion variable
!4 limitado por huevos con reabsorcion
!5 limitado por huevos con reabsorcion condicional

contains
  subroutine cambia_huevos_1(huevos_t, estimado)
    !limitado por huevos
    real*8 huevos_t,estimado
    if (estimado > huevos_t) then
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    return
  end subroutine cambia_huevos_1

  subroutine cambia_huevos_2(huevos_t,producci,estimado)
    real*8 huevos_t,producci,estimado
    if (estimado > huevos_t) then 
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    huevos_t = huevos_t + producci
    return
  end subroutine cambia_huevos_2

  subroutine cambia_huevos_3(huevos_t,producci,cambio_p,estimado)
    real*8 huevos_t,producci,cambio_p,estimado
    if (estimado > huevos_t) then 
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    huevos_t = huevos_t + producci
    producci = producci * cambio_p
    return
  end subroutine cambia_huevos_3

  subroutine cambia_huevos_4(huevos_t,producci,cambio_p,reabsorc,estimado)
    real*8 huevos_t,producci,cambio_p,reabsorc,estimado
    if (estimado > huevos_t) then 
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    huevos_t = huevos_t * reabsorc
    huevos_t = huevos_t + producci
    producci = producci * cambio_p
    return
  end subroutine cambia_huevos_4

  subroutine cambia_huevos_5(huevos_t,producci,cambio_p,reabsorc,umbral_r,estimado)
    real*8 huevos_t,producci,cambio_p,reabsorc,estimado,umbral_r
    if (estimado > huevos_t) then 
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    if (huevos_t > umbral_r) then
       huevos_t = huevos_t * reabsorc
    end if
    huevos_t = huevos_t + producci
    producci = producci * cambio_p
    return
  end subroutine cambia_huevos_5

  subroutine cambia_huevos_6(huevos_t,producci,cambio_p,reabsorc,capac_hu,estimado)
    real*8 huevos_t,producci,cambio_p,reabsorc,estimado,capac_hu
    if (estimado > huevos_t) then 
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    huevos_t = huevos_t * reabsorc
    huevos_t = huevos_t + producci
    if (huevos_t > capac_hu) then
       huevos_t = capac_hu
    end if
    producci = producci * cambio_p
    return
  end subroutine cambia_huevos_6

  subroutine cambia_huevos_7(huevos_t,producci,cambio_p,reabsorc,capac_hu,&
       umbral_r,estimado)
    real*8 huevos_t,producci,cambio_p,reabsorc,estimado,capac_hu,&
         umbral_r
    if (estimado > huevos_t) then 
       estimado = huevos_t
    end if
    huevos_t = huevos_t - estimado
    if (huevos_t > umbral_r) then
       huevos_t = huevos_t * reabsorc
    end if
    huevos_t = huevos_t + producci
    if (huevos_t > capac_hu) then
       huevos_t = capac_hu
    end if
    producci = producci * cambio_p
    return
  end subroutine cambia_huevos_7
  
  real*8 function modelobicho(dia,densidad,parasitados,params,&
       ndia,nden,npar,nparams,&
       modelo_huevos,resp_func)
    ! modelo de markov para cada bicho
    integer ndia, nden, npar, nparams, resp_func, modelo_huevos
    real*8 like, hp, he, nd, &
         dia(ndia),densidad(nden),parasitados(npar), &
         estimado, params(nparams), binomial, estimadoc
    real*8 huevos_t, producci, cambio_p, reabsorc, umbral_r, capac_hu, &
         infinity, cumhuevos
    parameter (infinity = 1.7976931348623157d300)
    like=0.
    hp=0.
    he=0.
    nd=0.
    cumhuevos = 0
    select case(modelo_huevos)
    case(1) 
       huevos_t = params(nparams)
    case(2)
       huevos_t = params(nparams-1)
       producci = params(nparams)
    case(3)
       huevos_t = params(nparams-2)
       producci = params(nparams-1)
       cambio_p = params(nparams)
    case(4)
       huevos_t = params(nparams-3)
       producci = params(nparams-2)
       cambio_p = params(nparams-1)
       reabsorc = params(nparams)
    case(5)
       huevos_t = params(nparams-4)
       producci = params(nparams-3)
       cambio_p = params(nparams-2)
       reabsorc = params(nparams-1)
       umbral_r = params(nparams)
    case(6)
       huevos_t = params(nparams-4)
       producci = params(nparams-3)
       cambio_p = params(nparams-2)
       reabsorc = params(nparams-1)
       capac_hu = params(nparams)
    case(7)
       huevos_t = params(nparams-5)
       producci = params(nparams-4)
       cambio_p = params(nparams-3)
       reabsorc = params(nparams-2)
       umbral_r = params(nparams-1)
       capac_hu = params(nparams)
    end select

    do i=1,ndia
       if (dia(i) == 0) exit
       if (densidad(i)==0) cycle
       cumhuevos = cumhuevos + densidad(i)
       select case(resp_func)
       case(1)
          estimado=respuesta_funcional_1(densidad(i),params(1))
       case(2)
          estimado=respuesta_funcional_2(densidad(i),params(2),params(1))
       case(3)
          estimado=respuesta_funcional_3(densidad(i),params(3),params(1),params(2))
       case(4)
          estimado=respuesta_funcional_4(densidad(i),params(3),params(1),params(2), cumhuevos)
       end select

       select case(modelo_huevos)
       case(1) 
          call cambia_huevos_1(huevos_t,estimado)
       case(2)
          call cambia_huevos_2(huevos_t,producci,estimado)
       case(3)
          call cambia_huevos_3(huevos_t,producci,cambio_p,estimado)
       case(4)
          call cambia_huevos_4(huevos_t,producci,cambio_p,reabsorc,estimado)
       case(5)
          call cambia_huevos_5(huevos_t,producci,cambio_p,reabsorc,umbral_r,estimado)
       case(6)
          call cambia_huevos_6(huevos_t,producci,cambio_p,reabsorc,capac_hu,estimado)
       case(7)
          call cambia_huevos_7(huevos_t,producci,cambio_p,reabsorc,capac_hu,umbral_r,estimado)

       end select
       estimadoc=estimado
       hp = hp + parasitados(i)
       he = he + estimado
       nd = nd + 1
       like = like + binomial(int(parasitados(i)), int(densidad(i)), estimado)
    end do

    modelobicho = like
    return
  end function modelobicho

  real*8 function modelobicho_v(dia,densidad,parasitados,params,&
       ndia,nden,npar,nparams,&
       modelo_huevos,resp_func,&
       salida)
    ! modelo de markov para cada bicho
    integer ndia, nden, npar, nparams, resp_func, modelo_huevos
    real*8 like, hp, he, nd, &
         dia(ndia),densidad(nden),parasitados(npar), &
         estimado, params(nparams), binomial, estimadoc
    real*8 huevos_t, producci, cambio_p, reabsorc, umbral_r, capac_hu, &
         infinity, cumhuevos
    real*8 salida(ndia)
    parameter (infinity = 1.7976931348623157d300)
    like=0.
    hp=0.
    he=0.
    nd=0.
    cumhuevos = 0
    select case(modelo_huevos)
    case(1) 
       huevos_t = params(nparams)
    case(2)
       huevos_t = params(nparams-1)
       producci = params(nparams)
    case(3)
       huevos_t = params(nparams-2)
       producci = params(nparams-1)
       cambio_p = params(nparams)
    case(4)
       huevos_t = params(nparams-3)
       producci = params(nparams-2)
       cambio_p = params(nparams-1)
       reabsorc = params(nparams)
    case(5)
       huevos_t = params(nparams-4)
       producci = params(nparams-3)
       cambio_p = params(nparams-2)
       reabsorc = params(nparams-1)
       umbral_r = params(nparams)
    case(6)
       huevos_t = params(nparams-4)
       producci = params(nparams-3)
       cambio_p = params(nparams-2)
       reabsorc = params(nparams-1)
       capac_hu = params(nparams)
    case(7)
       huevos_t = params(nparams-5)
       producci = params(nparams-4)
       cambio_p = params(nparams-3)
       reabsorc = params(nparams-2)
       umbral_r = params(nparams-1)
       capac_hu = params(nparams)
    end select

    do i=1,ndia
       if (dia(i) == 0) exit
       if (densidad(i)==0) cycle
       cumhuevos = cumhuevos + densidad(i)
       select case(resp_func)
       case(1)
          estimado=respuesta_funcional_1(densidad(i),params(1))
       case(2)
          estimado=respuesta_funcional_2(densidad(i),params(2),params(1))
       case(3)
          estimado=respuesta_funcional_3(densidad(i),params(3),params(1),params(2))
       case(4)
          estimado=respuesta_funcional_4(densidad(i),params(3),params(1),params(2), cumhuevos)
       end select

       select case(modelo_huevos)
       case(1) 
          call cambia_huevos_1(huevos_t,estimado)
       case(2)
          call cambia_huevos_2(huevos_t,producci,estimado)
       case(3)
          call cambia_huevos_3(huevos_t,producci,cambio_p,estimado)
       case(4)
          call cambia_huevos_4(huevos_t,producci,cambio_p,reabsorc,estimado)
       case(5)
          call cambia_huevos_5(huevos_t,producci,cambio_p,reabsorc,umbral_r,estimado)
       case(6)
          call cambia_huevos_6(huevos_t,producci,cambio_p,reabsorc,capac_hu,estimado)
       case(7)
          call cambia_huevos_7(huevos_t,producci,cambio_p,reabsorc,capac_hu,umbral_r,estimado)

       end select
       estimadoc = estimado
       hp = hp + parasitados(i)
       he = he + estimado
       nd = nd + 1
       salida(i) = estimado
    end do
    return
  end function modelobicho_v
end module modelos_bichos
