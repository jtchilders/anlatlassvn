module     p11_csbar_hepneg_dipoles
   use p11_csbar_hepneg_config, only: ki
   use p11_csbar_hepneg_color, only: gammaF, gammaA, &
       & CF, CA, numcs, KF, KA, &
       & T1T2, &
       & T1T6, &
       & T2T6
   use p11_csbar_hepneg_kinematics, only: num_legs, dotproduct, &
       lambda
   implicit none

   private :: ki, gammaF, gammaA, CF, CA, numcs, KF, KA, &
       & T1T2, &
       & T1T6, &
       & T2T6
   private :: num_legs, dotproduct, lambda


   real(ki), parameter :: pi = &
   & 3.1415926535897932384626433832795028841971693993751_ki

   private :: V_SING, GAMMA_F, GAMMA_A
   private :: I_ff, I_if, I_ii
contains
   function insertion_operator(mu_sq, vec, I, J)
      ! This function calculates the poles of the insertion operator
      ! for the specified process.
      ! See [Catani,Dittmaier,Seymour,Trocsanyi]
      ! The result does not include the factor of
      ! $-\alpha_s/(2\pi)(4\pi)^\epsilon/\Gamma(1-\epsilon)$
      implicit none
      real(ki), intent(in) :: mu_sq
      real(ki), dimension(num_legs,4), intent(in) :: vec
      integer, optional, intent(in) :: I, J
      complex(ki), dimension(numcs,numcs,2) :: insertion_operator

      if(present(I) .and. present(J)) then
              insertion_operator = &
                 &   I_ff(mu_sq, vec, I, J) &
                 & + I_if(mu_sq, vec, I, J) &
                 & + I_ii(mu_sq, vec, I, J)
      else
              insertion_operator = &
                 &   I_ff(mu_sq, vec) &
                 & + I_if(mu_sq, vec) &
                 & + I_ii(mu_sq, vec)
      end if
   end  function insertion_operator

   !=================================================================
   function I_ff(mu_sq, vec, I, J)
      ! This function calculates I_m as specified in (6.16) of
      ! [Catani,Dittmaier,Seymour,Trocsanyi]
      ! The result does not include the factor of
      ! $-\alpha_s/(2\pi)(4\pi)^\epsilon/\Gamma(1-\epsilon)$
      use p11_csbar_hepneg_model
      implicit none
      real(ki), intent(in) :: mu_sq
      real(ki), dimension(num_legs,4), intent(in) :: vec
      integer, optional, intent(in) :: I, J
      complex(ki), dimension(numcs,numcs,2) :: I_ff
      real(ki), dimension(2) :: term
      complex(ki), dimension(numcs,numcs) :: matrix
      real(ki) :: log_s, s_jk
      real(ki) :: Q_jk, rho, v_jk
      real(ki) :: m1_2, m2_2, mu1_2, mu2_2
      logical :: flag

      I_ff(:,:,:) = 0.0_ki
      
   end  function I_ff

   !=================================================================
   function I_if(mu_sq, vec, I, J)
      ! Implements equation (6.52) from the above paper for
      ! each initial state particle. Again we omit the prefactor
      ! $-\alpha_s/(2\pi)(4\pi)^\epsilon/\Gamma(1-\epsilon)$
      use p11_csbar_hepneg_model
      implicit none
      real(ki), intent(in) :: mu_sq
      real(ki), dimension(num_legs,4), intent(in) :: vec
      complex(ki), dimension(numcs,numcs,2) :: I_if
      real(ki), dimension(2) :: term
      complex(ki), dimension(numcs,numcs) :: matrix
      integer, optional, intent(in) :: I, J
      real(ki) :: log_s, s_ja
      logical :: flag

      I_if(:,:,:) = 0.0_ki

      s_ja = 2.0_ki * abs(dotproduct(vec(1,:), vec(6,:)))
      log_s = log(mu_sq / s_ja)
      !log_s = 0.0_ki

      ! T6.T1/T6^2
      if(present(I) .and. present(J)) then
         flag = (I .eq. 6) .and. (J .eq. 1)
      else
         flag = .true.
      end if

      if (flag) then
              matrix(:,:) = T1T6(:,:)/CA
              term = V_SING(s_ja, 0.0_ki, 0.0_ki)
              term(1) = term(1) + log_s * term(2)
              term = CA * term
              term(1) = term(1) + gammaA
              I_if(:,:,1) = I_if(:,:,1) + term(1) * matrix(:,:)
              I_if(:,:,2) = I_if(:,:,2) + term(2) * matrix(:,:)
      end if

      ! T1.T6/T1^2
      if(present(I) .and. present(J)) then
         flag = (I .eq. 1) .and. (J .eq. 6)
      else
         flag = .true.
      end if

      if (flag) then
              matrix(:,:) = T1T6(:,:)/CF

              term = V_SING(s_ja, 0.0_ki, 0.0_ki)
              term(1) = term(1) + log_s * term(2)
              term = CF * term
              term(1) = term(1) + gammaF
              I_if(:,:,1) = I_if(:,:,1) + term(1) * matrix(:,:)
              I_if(:,:,2) = I_if(:,:,2) + term(2) * matrix(:,:)
      end if

      s_ja = 2.0_ki * abs(dotproduct(vec(2,:), vec(6,:)))
      log_s = log(mu_sq / s_ja)
      !log_s = 0.0_ki

      ! T6.T2/T6^2
      if(present(I) .and. present(J)) then
         flag = (I .eq. 6) .and. (J .eq. 2)
      else
         flag = .true.
      end if

      if (flag) then
              matrix(:,:) = T2T6(:,:)/CA
              term = V_SING(s_ja, 0.0_ki, 0.0_ki)
              term(1) = term(1) + log_s * term(2)
              term = CA * term
              term(1) = term(1) + gammaA
              I_if(:,:,1) = I_if(:,:,1) + term(1) * matrix(:,:)
              I_if(:,:,2) = I_if(:,:,2) + term(2) * matrix(:,:)
      end if

      ! T2.T6/T2^2
      if(present(I) .and. present(J)) then
         flag = (I .eq. 2) .and. (J .eq. 6)
      else
         flag = .true.
      end if

      if (flag) then
              matrix(:,:) = T2T6(:,:)/CF

              term = V_SING(s_ja, 0.0_ki, 0.0_ki)
              term(1) = term(1) + log_s * term(2)
              term = CF * term
              term(1) = term(1) + gammaF
              I_if(:,:,1) = I_if(:,:,1) + term(1) * matrix(:,:)
              I_if(:,:,2) = I_if(:,:,2) + term(2) * matrix(:,:)
      end if
   end  function I_if

   !=================================================================
   function I_ii(mu_sq, vec, I, J)
      ! This function covers those terms from eq. (6.66)
      ! which are not part of the other two fuctions
      ! The result omits the factor
      ! $-\alpha_s/(2\pi)(4\pi)^\epsilon/\Gamma(1-\epsilon)$
      implicit none
      real(ki), intent(in) :: mu_sq
      real(ki), dimension(num_legs,4), intent(in) :: vec
      integer, optional, intent(in) :: I, J
      complex(ki), dimension(numcs,numcs,2) :: I_ii
      real(ki), dimension(2) :: term
      complex(ki), dimension(numcs,numcs) :: matrix
      real(ki) :: log_s, s_ab
      logical :: flag

      I_ii(:,:,:) = 0.0_ki
      
      s_ab = 2.0_ki * abs(dotproduct(vec(1,:), vec(2,:)))
      log_s = log(mu_sq / s_ab)
      !log_s = 0.0_ki

      ! T1.T2/T1^2
      if(present(I) .and. present(J)) then
         flag = (I .eq. 1) .and. (J .eq. 2)
      else
         flag = .true.
      end if

      if (flag) then
              matrix(:,:) = T1T2(:,:)/CF
              term(1) = gammaF
              term(2) = CF
              term(1) = term(1) + log_s * term(2)

              I_ii(:,:,1) = I_ii(:,:,1) + term(1) * matrix(:,:)
              I_ii(:,:,2) = I_ii(:,:,2) + term(2) * matrix(:,:)
      end if

      ! T2.T1/T2^2
      if(present(I) .and. present(J)) then
         flag = (I .eq. 2) .and. (J .eq. 1)
      else
         flag = .true.
      end if

      if (flag) then
              matrix(:,:) = T1T2(:,:)/CF
              term(1) = gammaF
              term(2) = CF
              term(1) = term(1) + log_s * term(2)

              I_ii(:,:,1) = I_ii(:,:,1) + term(1) * matrix(:,:)
              I_ii(:,:,2) = I_ii(:,:,2) + term(2) * matrix(:,:)
      end if

   end  function I_ii

   function V_SING(s_jk, m_j, m_k)
      implicit none
      real(ki), intent(in) :: s_jk, m_j, m_k
      real(ki), dimension(2) :: V_SING

      real(ki) :: Q_jk, v_jk, rho, mu1_2, mu2_2, m1_2, m2_2

      m1_2  = m_j*m_j
      m2_2  = m_k*m_k

      if (m_j .gt. 0.0_ki) then
              if(m_k .gt. 0.0_ki) then
                      Q_jk  = s_jk + m1_2 + m2_2
                      mu1_2 = m1_2 / Q_jk
                      mu2_2 = m2_2 / Q_jk
                      v_jk  = sqrt(lambda(1.0_ki, mu1_2, mu2_2)) &
                            &      / (1.0_ki - mu1_2 - mu2_2)
                      rho   = sqrt((1.0_ki - v_jk)/(1.0_ki + v_jk))
                      ! (6.20), case 1, both massive
                      V_SING(1) = log(rho)/v_jk
                      V_SING(2) = 0.0_ki
              else
                      ! (6.20), case 2, only m_j massive
                      V_SING(1) = 0.5_ki * log(m1_2/s_jk)
                      V_SING(2) = 0.5_ki
              end if
      else
              if(m_k .gt. 0.0_ki) then
                      ! (6.20), case 2, only m_k massive
                      V_SING(1) = 0.5_ki * log(m2_2/s_jk)
                      V_SING(2) = 0.5_ki
              else
                      ! (6.20), case 3, both massless
                      V_SING(1) = 0.0_ki
                      V_SING(2) = 1.0_ki
              end if
      end if
   end  function V_SING

   function GAMMA_A()
      ! Equation (6.27)
      implicit none
      real(ki), dimension(2) :: GAMMA_A

      GAMMA_A(1) = gammaA
      GAMMA_A(2) = 0.0_ki
   end  function GAMMA_A

   function GAMMA_F(mq)
      implicit none
      real(ki), intent(in) :: mq
      real(ki), dimension(2) :: GAMMA_F

      if (mq .gt. 0.0_ki) then
              ! (6.29)
              GAMMA_F(1) = CF
              GAMMA_F(2) = 0.0_ki
      else
              ! (6.28)
              GAMMA_F(1) = gammaF
              GAMMA_F(2) = 0.0_ki
      end if
   end  function GAMMA_F
end module p11_csbar_hepneg_dipoles
