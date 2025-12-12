program simple
  use ziggurat
  use globals
  use md_module
  use md_integrator
  implicit none

  integer :: seed
  logical :: es
  integer :: step, step_eq, step_prod
  real(kind=8) :: K, Tinst, Etotal, P, Pmed, Pmed2, conteoVirial, var_presion
  real(kind=8), allocatable :: f_old(:,:)

!=========================================================  
  ! inicializacion
!=========================================================
  call init()
  allocate(f_old(3,N))
  Pmed =  0
  Pmed2 = 0
  conteoVirial = 0
!=========================================================
! Semilla
!=========================================================
    inquire(file='seed.dat', exist=es)
    if (es) then
      open(10, file='seed.dat', status='old')
      read(10,*) seed
      close(10)
      print *, " * Leyendo semilla de archivo seed.dat"
    else
      seed = 1234567
    end if
    call zigset(seed)
!=========================================================
! inicio posicion y velocidades. Calulco la fuerza
!=========================================================
    call init_positions()
    call apply_pbc_positions()
    call init_velocities(Ttarget)

    call force(f, E_pot, W_inst)
    
    !abro archivos outputs
    open(unit=50, file='output/energia.dat', status='unknown')
    open(unit=20,file="output/posicionT.xyz",status="unknown")
    open(unit=40,file="output/presion.dat",status="unknown")
    open(unit=30,file="output/temperatura.dat",status="unknown")
    open(unit=10,file="output/fluctP.dat",status="unknown")
    open(unit=60,file="output/varP.dat",status="unknown")
    
    !preparar muestreo g(r)
    call init_gr(ngr_bins, rmax)  
    
!=========================================================
! Minimización de la energía
!=========================================================

  print *, "Iniciando minimización. E_pot inicial:", E_pot

  do step = 1, tMinimizacion
    r(:,:) = r(:,:) + 0.5d0 * f(:,:) * deltaTMinimizacion * deltaTMinimizacion
    call apply_pbc_positions()
    call force(f, E_pot, W_inst)
    ! Llamada a wall93 para calcular fuerzas fluido-pared
    call wall93(2, r, f, a_wall, sigma_wall, L, N, a_type)

    if (mod(step, muestreoMinimizacion)==0) then
             Etotal = E_pot
             write(50, *) Etotal
    end if
  end do
  
!=========================================================
! Termalización
!=========================================================

    do step_eq = 1, 100000
        f_old(:,:) = f(:,:)
        call verlet_positions(deltaT)
        call apply_pbc_positions()
        call force(f, E_pot, W_inst)
        call add_langevin_forces(f, v, Ttarget, gama, deltaT)
        ! Llamada a wall93 para calcular fuerzas fluido-pared
        call wall93(2, r, f, a_wall, sigma_wall, L, N, a_type)
        call verlet_velocities(deltaT, f_old)

        if (mod(step_eq, muestreo) == 0) then
            call measurements(K, Tinst, P, Etotal)  
            write(50, *) Etotal     
            write(40, *) P
            write(30,*) Tinst
            call write_xyz(20, step_prod, Etotal)
        end if
        
    end do
    print *, "Equilibración finalizada."

!=========================================================
! Producción (y muestreo)
!=========================================================
    do step_prod = 1, pasosT - tMinimizacion
        f_old(:,:) = f(:,:)
        call verlet_positions(deltaT)
        call apply_pbc_positions()
        call force(f, E_pot, W_inst)
        call add_langevin_forces(f, v, Ttarget, gama, deltaT)
        call add_pumping(f, pump, N)
        ! Llamada a wall93 para calcular fuerzas fluido-pared
        call wall93(2, r, f, a_wall, sigma_wall, L, N, a_type)
        call verlet_velocities(deltaT, f_old)
        
        if (mod(step_prod, muestreo) == 0) then
            call measurements(K, Tinst, P, Etotal)  
            Pmed =  Pmed + P
            Pmed2 = Pmed2 + P*P
            conteoVirial= conteoVirial + W_inst
            write(50, *) Etotal     
            write(40, *) P
            write(30,*) Tinst
            call write_xyz(20, step_prod, Etotal)
        end if
        if (mod(step_prod,1000) == 0) call accumulate_gr() 

        ! Imprimir porcentaje de avance
        if (mod(step_prod,int((pasosT - tMinimizacion)/100)) == 0) then 
            print '(A,F6.2,A)', " Progreso: ", 100.0d0*step_prod/(pasosT - tMinimizacion), "%"
        end if  
     
     end do
!===============================
! posprocesamiento
!===============================

    Pmed   =  Pmed / (REAL(pasosT - tMinimizacion)/REAL(MUESTREO))
    Pmed2  = Pmed2 / (REAL(pasosT - tMinimizacion)/REAL(MUESTREO))
    conteoVirial = conteoVirial / (REAL(pasosT - tMinimizacion)/REAL(MUESTREO))
    var_presion = sqrt(Pmed2 - (Pmed*Pmed))
    write(10, *) Pmed2, Pmed
    write(60,*) var_presion
    
    close(20)
    close(40)
    close(50)
    close(30)
    close(10)
    close(60)
    !Normalizar y escribir g(r)
    call normalize_and_write_gr("gr.dat")
    print *, "g(r) escrito en gr.dat"

!===============================
! Guardar nueva semilla
!===============================
    open(10,file='seed.dat',status='unknown')
    seed = shr3()
    write(10,*) seed
    close(10)

    print *, "MD (Guia 5) finalizado."

end program simple

