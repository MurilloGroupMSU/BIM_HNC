! This solves HNC equations for a binary mixtures (BIM)

! Length units are ion-sphere radii
! Inputs are gamma=e^2/(a T), Z1, Z2, x

program BIM_HNC

    implicit none

    integer, parameter :: num_bins=300, num_iter =100
    real, parameter :: r_max = 10.0, Four_Pi = 12.5663706144
    real, parameter :: Pi = 3.141592653589793, Th_O_F_Pi=0.238732415
    real, parameter :: Two_Pi_Sq = 19.7392088022, ln_2 = 0.693147182
    real, dimension(num_bins) :: u_s_11_r, u_s_22_r, u_s_12_r, u_s_21_r
    real, dimension(num_bins) :: u_l_11_k, u_l_22_k, u_l_12_k, u_l_21_k
    real, dimension(num_bins) :: c_11_k, c_22_k, c_12_k, c_21_k
    real, dimension(num_bins) :: c_s_11_k, c_s_22_k, c_s_12_k, c_s_21_k
    real, dimension(num_bins) :: c_s_11_r, c_s_22_r, c_s_12_r, c_s_21_r
    real, dimension(num_bins) :: N_s_11_k, N_s_22_k, N_s_12_k, N_s_21_k
    real, dimension(num_bins) :: N_s_11_r, N_s_22_r, N_s_12_r, N_s_21_r
    real, dimension(num_bins) :: g_11, g_22, g_21, g_12
    real, dimension(num_bins) :: b_r_11, b_r_22, b_r_12,b_r_21
    real :: gamma, n_1, n_2, n_e
    real :: gamma_11, gamma_22, gamma_12, gamma_21
    real :: r, del_r, k, del_k, alpha, determ
    integer :: bin, it_loop, bin_r, bin_k
    real :: x, charge_1,charge_2, Bridge_status
    ! Yukawa bridge function
!real :: C11, C22,C12
    real :: b_0, b_1, c_1, c_2, c_3, kappa, step
    real:: a11,a22,a12,a21
    real :: x1,x2

  !  open(unit=13, file='dcfs.out', status='unknown')
    open(unit=14, file='rdfs.out', status='unknown')


 ! Inputs parameters using bim.in
 
     ! x=ni/n species isotopic ratio. ni is i-species density
     ! n is the total density.
   open(unit=1, status='unknown', file='bim.in')
  
   read(1,*) gamma
   read(1,*) x
   read(1,*) charge_1
   read(1,*) charge_2
   read(1,*) Bridge_status

   ! gamma=10.271    ! e^2/aT where a is the total Wigner-Setiz radius.
   ! x=0.05          !  
   ! Z1=1.           ! Species 1 charge 
    !charge_2=5.           ! Species 2 charge 
    !Brige_off
    
    
    
    n_1 = (1.-x)*Th_O_F_Pi
    n_2 =x*Th_O_F_Pi
    
! gamma_ij= Zi*Zj/aT
    
    gamma_11 = charge_1**2.*gamma    
    gamma_22 = charge_2**2.*gamma
    gamma_12 = charge_1*charge_2*gamma
    gamma_21 = charge_1*charge_2*gamma
    
    kappa =  0. 
    alpha =1. 

! Steps in physical and Fourier spaces

    del_r = r_max/num_bins
    del_k = Pi/r_max


    do bin = 1, num_bins

        r = bin*del_r
        k = bin*del_k
     
           
    ! Bridge function. The coeffiencts b_1 and b_1 are calculated using Iyetomi, et al., PRA 46, 1051 (1992)
       
  
       b_r_11(bin) =-.0464*gamma_11**(1.336)*exp(-b_1(gamma_11)/b_0(gamma_11)*(r/1.25)**2.)
       b_r_22(bin) =-.0464*gamma_22**(1.336)*exp(-b_1(gamma_22)/b_0(gamma_22)*(r/1.25)**2.)
       b_r_21(bin) =-.0464*gamma_21**(1.336)*exp(-b_1(gamma_21)/b_0(gamma_21)*(r/1.25)**2.)
       b_r_12(bin) =-.0464*gamma_12**(1.336)*exp(-b_1(gamma_12)/b_0(gamma_12)*(r/1.25)**2.)
        
    
   !effective short-ranged potential

        u_s_11_r(bin) = gamma_11*(exp(-(alpha+kappa)*r))/r
        u_s_12_r(bin) = gamma_12*(exp(-(alpha+kappa)*r))/r
        u_s_21_r(bin) = gamma_21*(exp(-(alpha+kappa)*r))/r
        u_s_22_r(bin) = gamma_22*(exp(-(alpha+kappa)*r))/r


    ! effective long-ranged potential

       u_l_11_k(bin) = Four_Pi*gamma_11*alpha**2/(k**2*(k**2+alpha**2))
       u_l_12_k(bin) = Four_Pi*gamma_12*alpha**2/(k**2*(k**2+alpha**2))
       u_l_21_k(bin) = Four_Pi*gamma_21*alpha**2/(k**2*(k**2+alpha**2))
       u_l_22_k(bin) = Four_Pi*gamma_22*alpha**2/(k**2*(k**2+alpha**2))


    ! Intial guess for the direct correlation functions
    
       c_11_k(bin) = -four_pi*gamma_11/(k**2+kappa**2) ! 
       c_12_k(bin) = -four_pi*gamma_12/(k**2+kappa**2) !
       c_21_k(bin) = -four_pi*gamma_21/(k**2+kappa**2) !
       c_22_k(bin) = -four_pi*gamma_22/(k**2+kappa**2) !
       
       c_s_11_k(bin) = c_11_k(bin) + u_l_11_k(bin)
       c_s_12_k(bin) = c_12_k(bin) + u_l_12_k(bin)
       c_s_21_k(bin) = c_21_k(bin) + u_l_21_k(bin)
       c_s_22_k(bin) = c_22_k(bin) + u_l_22_k(bin)

       N_s_11_k(bin) = 0.0
       N_s_22_k(bin) = 0.0
       N_s_12_k(bin) = 0.0
       N_s_21_k(bin) = 0.0

       N_s_11_r(bin) = 0.0
       N_s_22_r(bin) = 0.0
       N_s_12_r(bin) = 0.0
       N_s_21_r(bin) = 0.0

    enddo

! Main Iteration Loop

    do it_loop = 1, num_iter



        ! Find N_s(k)
        do bin = 1, num_bins

        determ = (1.0 - n_1*(c_s_11_k(bin)-u_l_11_k(bin)))* &
        & (1.0 - n_2*(c_s_22_k(bin)-u_l_22_k(bin))) - n_1*n_2* &
        & (c_s_12_k(bin)-u_l_12_k(bin))*(c_s_21_k(bin)-u_l_21_k(bin))


        N_s_11_k(bin) = ( (c_s_11_k(bin)-u_l_11_k(bin)))/&
        & determ - c_s_11_k(bin)

        N_s_11_k(bin) = ( (c_s_11_k(bin)-u_l_11_k(bin))*(1.0 - &
        & n_2*(c_s_22_k(bin)-u_l_22_k(bin))) + n_2*&
        & (c_s_12_k(bin)-u_l_12_k(bin))*(c_s_21_k(bin)-u_l_21_k(bin)) )/&
        & determ - c_s_11_k(bin)

        N_s_22_k(bin) = ( (c_s_22_k(bin)-u_l_22_k(bin))*(1.0 - &
        & n_1*(c_s_11_k(bin)-u_l_11_k(bin))) + n_1*&
        & (c_s_12_k(bin)-u_l_12_k(bin))*(c_s_12_k(bin)-u_l_12_k(bin)) )/&
        & determ - c_s_22_k(bin)

        N_s_12_k(bin) = ( (c_s_12_k(bin)-u_l_12_k(bin))*(1.0 - &
        & n_2*(c_s_22_k(bin)-u_l_22_k(bin))) + n_2* &
        & (c_s_12_k(bin)-u_l_12_k(bin))*(c_s_22_k(bin)-u_l_22_k(bin)) )/&
        & determ - c_s_12_k(bin)

        N_s_21_k(bin) = ( (c_s_21_k(bin)-u_l_21_k(bin))*(1.0 - &
        & n_1*(c_s_11_k(bin)-u_l_11_k(bin))) + n_1* &
        & (c_s_21_k(bin)-u_l_21_k(bin))*(c_s_11_k(bin)-u_l_11_k(bin)) )/&
        & determ - c_s_21_k(bin)

       enddo
       
    ! Find N_s(r)
    
        do bin_r = 1, num_bins

            r = bin_r*del_r
            N_s_11_r(bin_r) = 0.0
            N_s_22_r(bin_r) = 0.0
            N_s_21_r(bin_r) = 0.0
            N_s_12_r(bin_r) = 0.0

                    do bin_k = 1, num_bins

                        k = bin_k*del_k

                        N_s_11_r(bin_r) = N_s_11_r(bin_r) + k*sin(k*r)*N_s_11_k(bin_k)
                        N_s_22_r(bin_r) = N_s_22_r(bin_r) + k*sin(k*r)*N_s_22_k(bin_k)
                        N_s_21_r(bin_r) = N_s_21_r(bin_r) + k*sin(k*r)*N_s_21_k(bin_k)
                        N_s_12_r(bin_r) = N_s_12_r(bin_r) + k*sin(k*r)*N_s_12_k(bin_k)

                    end do

            N_s_11_r(bin_r) = N_s_11_r(bin_r)*del_k/(Two_Pi_Sq*r)
            N_s_22_r(bin_r) = N_s_22_r(bin_r)*del_k/(Two_Pi_Sq*r)
            N_s_21_r(bin_r) = N_s_21_r(bin_r)*del_k/(Two_Pi_Sq*r)
            N_s_12_r(bin_r) = N_s_12_r(bin_r)*del_k/(Two_Pi_Sq*r)

     end do
     
    ! Find g(r) and C_s(r)
    
     do bin = 1, num_bins

        g_11(bin) = exp(-u_s_11_r(bin) + N_s_11_r(bin)+ Bridge_status*b_r_11(bin))
        g_22(bin) = exp(-u_s_22_r(bin) + N_s_22_r(bin)+ Bridge_status*b_r_22(bin))
        g_21(bin) = exp(-u_s_21_r(bin) + N_s_21_r(bin)+ Bridge_status*b_r_21(bin))
        g_12(bin) = exp(-u_s_12_r(bin) + N_s_12_r(bin)+ Bridge_status*b_r_12(bin))

        c_s_11_r(bin) = g_11(bin) - 1.0 - N_s_11_r(bin)
        c_s_22_r(bin) = g_22(bin) - 1.0 - N_s_22_r(bin)
        c_s_21_r(bin) = g_21(bin) - 1.0 - N_s_21_r(bin)
        c_s_12_r(bin) = g_12(bin) - 1.0 - N_s_12_r(bin)

     end do


        ! Find c_s(k)

        do bin_k = 1, num_bins

            k = bin_k*del_k
            c_s_11_k(bin_k) = 0.0
            c_s_22_k(bin_k) = 0.0
            c_s_21_k(bin_k) = 0.0
            c_s_12_k(bin_k) = 0.0

            do bin_r = 1, num_bins

            r = bin_r*del_r

            c_s_11_k(bin_k) = c_s_11_k(bin_k) + r*sin(k*r)*c_s_11_r(bin_r)
            c_s_22_k(bin_k) = c_s_22_k(bin_k) + r*sin(k*r)*c_s_22_r(bin_r)
            c_s_21_k(bin_k) = c_s_21_k(bin_k) + r*sin(k*r)*c_s_21_r(bin_r)
            c_s_12_k(bin_k) = c_s_12_k(bin_k) + r*sin(k*r)*c_s_12_r(bin_r)

            end do

            c_s_11_k(bin_k) = c_s_11_k(bin_k)*Four_Pi*del_r/k
            c_s_22_k(bin_k) = c_s_22_k(bin_k)*Four_Pi*del_r/k
            c_s_21_k(bin_k) = c_s_21_k(bin_k)*Four_Pi*del_r/k
            c_s_12_k(bin_k) = c_s_12_k(bin_k)*Four_Pi*del_r/k
            
       !     write(13,*)  k, c_11_k(bin_k), c_12_k(bin_k), c_22_k(bin_k)

        end do

    end do

! Write final results radial pair distribution functions

        do bin = 1, num_bins

            k = bin*del_k
            r = bin*del_r
          
            write(14,*) r,  g_11(bin),g_22(bin),g_12(bin)
         enddo

end program BIM_HNC



    ! Bridge function parameters (Iyetomi, et al., PRA 46, 1051 (1992))
    ! Do not use for gamma<5.0


    FUNCTION b_0(gamma)
        real gamma
        real b_0

        b_0 = 0.258 - 0.0612*log(gamma) + 0.0123*(log(gamma))**2 - 1./gamma
        return

    End Function b_0

    FUNCTION b_1(gamma)
    real gamma
    real b_1

    b_1 = 0.0269 + 0.0318*log(gamma) + 0.00814*(log(gamma))**2
    return

    End Function b_1

    FUNCTION c_1(gamma)
    real gamma
    real c_1

    c_1 = 0.498 - 0.28*log(gamma) + 0.0294*(log(gamma))**2
    return

    End Function c_1

    FUNCTION c_2(gamma)
    real gamma
    real c_2

    c_2 = -0.412 + 0.219*log(gamma) - 0.0251*(log(gamma))**2
    return

    End Function c_2
    FUNCTION c_3(gamma)
    real gamma
    real c_3

    c_3 = 0.0988 - 0.0534*log(gamma) + 0.00682*(log(gamma))**2
    return

    End Function c_3


