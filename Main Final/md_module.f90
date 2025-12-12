module md_module
  use globals
  implicit none
contains

  !==============================================
  !Rutina de inicio
  !==============================================
  
  subroutine init()
    use ziggurat
    implicit none

    !Lectura de parámetros
    N        = Lectura_de_parametros('input/N')
    N_plus   = Lectura_de_parametros('input/N_plus')
    sigma    = Lectura_de_parametros('input/sigma')
    pump     = Lectura_de_parametros('input/pump')
    epsilon  = 1.0d0
    
    deltaTMinimizacion = Lectura_de_parametros('input/deltaTMinimizacion')
    deltaT             = Lectura_de_parametros('input/deltaT')
    pasosT             = int(Lectura_de_parametros('input/pasosT'))
    tMinimizacion      = int(Lectura_de_parametros('input/tMinimizacion'))
    muestreo           = int(Lectura_de_parametros('input/muestreo'))
    muestreoMinimizacion = int(Lectura_de_parametros('input/muestreoMinimizacion'))
    
    ! Inicialización de parámetros de interacción fluido-pared
    a_wall(1) = 1   ! Coeficiente para la pared superior
    a_wall(2) = 1   ! Coeficiente para la pared inferior
    sigma_wall = 1  ! Parámetro de interacción partícula-pared
    a_type = 1      ! Tipo de partícula (único tipo en este caso)

    Ttarget = Lectura_de_parametros('input/Temp')
    gama    = 0.5d0

    !caja
    L = Lectura_de_parametros('input/L')

    ! g(r) parámetros para función de distribución radial
    ngr_bins = 100
    rmax     = 0.5d0 * L

    ! allocate arrays
    if (.not.allocated(r)) allocate(r(3,N+N_plus))
    if (.not.allocated(v)) allocate(v(3,N))
    if (.not.allocated(f)) allocate(f(3,N+N_plus))
    
    E_pot = 0d0
    W_inst = 0d0
    f(:,:) = 0d0
    v(:,:) = 0d0

  end subroutine init
  
  !==========================================================
  ! Fuerza, energía potencial virial (para calculo de presión)
  !==========================================================

subroutine force(f_out, en, virial_local)
  use globals
  implicit none
  real(kind=8), intent(out) :: f_out(3,N+N_plus)
  real(kind=8), intent(out) :: en
  real(kind=8), intent(out) :: virial_local
  real(kind=8) :: deltar(3)
  real(kind=8) :: rij2, rij
  real(kind=8) :: sr2, sr6, sr12
  real(kind=8) :: f_mod, V_LJ
  integer :: i, j, k
  
  f_out(:,:) = 0.0d0
  en = 0.0d0
  virial_local = 0.0d0
  
  rc = 2.5d0 * sigma
  rc2 = rc*rc
  
  !para igualar a cero el potencial de corte
  V_shift = 4d0*epsilon*((sigma/rc)**12 - (sigma/rc)**6)
  
  do j = 1, N-1
    do i = j+1, N+N_plus
      deltar = r(:,i) - r(:,j)
      
      !Mìnima imagen
      do k = 1, 3
        if (deltar(k) > L/2.0d0) then
          deltar(k) = deltar(k) - L     ! Si la distancia es mayor que L/2, usar imagen periódica
        else if (deltar(k) < -L/2.0d0) then
          deltar(k) = deltar(k) + L
        end if
      end do
      
      rij2 = deltar(1)*deltar(1) + deltar(2)*deltar(2) + deltar(3)*deltar(3)
      
      ! Verificar que no estamos en la misma posición
      if (rij2 < 1d-12) then
        print *, "ADVERTENCIA: Partículas", i, "y", j, "en la misma posición"
        cycle
      end if
      
      ! Verificar cutoff
      if (rij2 > rc2) cycle
      
               ! Cálculo de energía y fuerza de Lennard-Jones
               rij = sqrt(rij2)
               sr2 = (sigma/rij)**2
               sr6 = sr2 * sr2 * sr2
               sr12 = sr6 * sr6
      
              ! Energía potencial
              V_LJ = 4.0d0 * epsilon * (sr12 - sr6)
              en = en + (V_LJ - V_shift)
       
              ! Fuerza: F = -dV/dr
              f_mod = 24.0d0 * epsilon * (2.0d0 * sr12 - sr6) / rij2
      
              !Aplico fuerzas
              f_out(:,i) = f_out(:,i) + f_mod * deltar
              f_out(:,j) = f_out(:,j) - f_mod * deltar
      
              ! Contribución al virial
              virial_local = virial_local + f_mod * rij2
    end do
  end do
  
  E_pot = en
  W_inst = virial_local
end subroutine force



  !==========================================================
  !fuerza de langevin, termostato
  !==========================================================

subroutine add_langevin_forces(f_tot, vel, Tset, gama_L, dt_local)
    use globals
    use ziggurat
    implicit none

    real(kind=8), intent(inout) :: f_tot(3,N)   ! fuerzas ya calculadas
    real(kind=8), intent(in)    :: vel(3,N)     ! velocidades actuales
    real(kind=8), intent(in)    :: Tset         ! temperatura objetivo
    real(kind=8), intent(in)    :: gama_L       ! coef de friccion
    real(kind=8), intent(in)    :: dt_local     ! delta t

    real(kind=8) :: ruido_amp
    integer :: i

    ! ---- amplitud del ruido ----
    ruido_amp = sqrt( 2d0 * gama_L * Tset / dt_local )

    do i = 1, N

        ! fuerza de friccion
        f_tot(:,i) = f_tot(:,i) - gama_L * vel(:,i)

        ! ruido gaussiano independiente por coordenada
        ! f_tot(1,i) = f_tot(1,i) + ruido_amp * rnor()
        f_tot(2,i) = f_tot(2,i) + ruido_amp * rnor()
        f_tot(3,i) = f_tot(3,i) + ruido_amp * rnor()

    end do

end subroutine add_langevin_forces

  !==========================================================
  !función para la lecturra de parámetros
  !==========================================================

 FUNCTION Lectura_de_parametros(nombre_archivo) RESULT (dato)
       CHARACTER(LEN=*), INTENT(IN) :: nombre_archivo
       real (kind=8) :: dato
    
       open(unit=10,file=nombre_archivo,status='old')
           read(10,*) dato
       close(10)
 END FUNCTION
 
  !==========================================================
  !escribe la configuracion de posiciones
  !==========================================================
 
 subroutine write_xyz(unit, step, Etotal)
    use globals
    implicit none
    integer, intent(in) :: unit
    integer, intent(in) :: step
    real(kind=8), intent(in) :: Etotal

    integer :: j, k
    real(kind=8) :: r_out(3)

    write(unit,*) N+N_plus
    write(unit,*) "Step=", step, " Etot=", Etotal

    do j = 1, N+N_plus
        !aplicar PBC a cada coordenada
        do k = 1, 3
            r_out(k) = r(k,j)
            if (r_out(k) >= L) r_out(k) = r_out(k) - L
            if (r_out(k) <  0d0) r_out(k) = r_out(k) + L
        end do

        write(unit,*) "A", r_out(1), r_out(2), r_out(3)
    end do
end subroutine write_xyz

  !==========================================================
  !Agrega las fuerzas de bombeo
  !==========================================================
  subroutine add_pumping(f_in, pump_force, n_particles)
      implicit none

      integer, intent(in)         :: n_particles
      real(kind=8), intent(inout) :: f_in(3,n_particles)   ! fuerzas ya calculadas
      real(kind=8), intent(in)    :: pump_force            ! fuerza de bombeo
      integer :: i

      do i = 1, n_particles
        f_in(1,i) = f_in(1,i) + pump_force
      end do

  end subroutine add_pumping

  !===========================================================
  ! Wall 9-3
  !===========================================================
  
 subroutine wall93(inter_type, r0, force, a_wall, sigma_wall, z_space_wall, n_mon_tot, a_type)
  ! Computes fluid-wall interactions:
  ! inter_type: Tipo de interacción (1, 2, etc.)
  ! r0: Posiciones de las partículas (3 x n_mon_tot)
  ! force: Fuerzas en las partículas (3 x n_mon_tot)
  ! a_wall: Coeficientes de interacción partícula-pared (2 x tipos)
  ! sigma_wall: Parámetro de interacción partícula-pared (tipos)
  ! z_space_wall: Altura del canal
  ! n_mon_tot: Número total de partículas
  ! a_type: Tipo de cada partícula (1 x n_mon_tot)

  implicit none
  integer, intent(in) :: inter_type, n_mon_tot
  real(kind=8), intent(in) :: r0(3, n_mon_tot), z_space_wall
  real(kind=8), intent(in) :: a_wall(2), sigma_wall
  integer, intent(in) :: a_type
  real(kind=8), intent(inout) :: force(3, n_mon_tot)

  real(kind=8) :: inv_z, r_dummy, v_fluid_wall
  integer :: i_part, i_type

  v_fluid_wall = 0.0d0
  
  select case(inter_type)
    case(2) ! (1/z)^9 - (1/z)^3 interaction
      do i_part = 1, n_mon_tot

        ! *** Bottom wall interaction
        i_type = a_type
        inv_z = 1.0d0 / r0(3, i_part)

        ! Limitamos la cercania de la partícula a la pared 
        if (inv_z > 20) inv_z = 20 

        r_dummy = sigma_wall * inv_z
        v_fluid_wall = v_fluid_wall + abs(a_wall(2)) * r_dummy**9 - a_wall(2) * r_dummy**3
        force(3, i_part) = force(3, i_part) + 9.0d0 * abs(a_wall(2)) * (sigma_wall)**9 * (inv_z)**10
        force(3, i_part) = force(3, i_part) - 3.0d0 * a_wall(2) * (sigma_wall)**3 * (inv_z)**4

        ! *** Top wall interaction
        inv_z = 1.0d0 / (z_space_wall - r0(3, i_part))
        
        ! Limitamos la cercania de la partícula a la pared 
        if (inv_z > 20) inv_z = 20 

        r_dummy = sigma_wall * inv_z
        v_fluid_wall = v_fluid_wall + abs(a_wall(1)) * r_dummy**9 - a_wall(1) * r_dummy**3
        force(3, i_part) = force(3, i_part) - 9.0d0 * abs(a_wall(1)) * (sigma_wall)**9 * (inv_z)**10
        force(3, i_part) = force(3, i_part) + 3.0d0 * a_wall(1) * (sigma_wall)**3 * (inv_z)**4

      end do

    case default
      print *, "Error: inter_type debe ser 1, 2, 3 o 4."
      !stop
  end select
end subroutine wall93

end module md_module



