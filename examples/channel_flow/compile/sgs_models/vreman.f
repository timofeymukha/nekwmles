c-----------------------------------------------------------------------
c> Compute eddy viscosity using Vreman model
c>  Here is the key difference:
c>  - no wall damping is required
c>  - algebraic, no dynamic procedure
c>
c> Refer: Vreman, A. W., Physics of Fluids (2004): https://doi.org/10.1063/1.1785131.
c> @callgraph @callergraph
      subroutine eddy_visc(e)
      implicit none

      include 'SIZE'
      include 'INPUT'  ! param
      include 'SOLN'  ! vx
      include 'SGS'  ! ediff, sij, snrm, dg2, dg2_max

      integer e
      real C_vreman, delta_sq, alpha_sq
      integer i, j, k, ntot

      real alpha(ldim, ldim), beta(ldim, ldim)
      common /vreman/ alpha, beta
      real wmles_sgs_c0

      ntot = nx1*ny1*nz1

c------need to be by element ->
      call comp_gije(sij, vx(1,1,1,e), vy(1,1,1,e), vz(1,1,1,e), e)
      ! Set gradient to match boundary condition

      wmles_sgs_c0 = 0.16
      C_vreman = 2.5 * wmles_sgs_c0**2

      do i=1, ntot
        delta_sq = dg2_max(e)

        do k=1,ldim
          do j=1,ldim
            alpha(j, k) = sij(i, j, k)
          end do
        end do

        beta = delta_sq * matmul(transpose(alpha), alpha)

        !> Note \f$ B_\beta \f$ is being stored in snrm here
        snrm(i, e) = (
     &      beta(1,1) * beta(2,2)
     &    - beta(1,2)**2  ! check
     &    + beta(1,1) * beta(3,3)
     &    - beta(1,3)**2  ! check
     &    + beta(2,2) * beta(3,3)
     &    - beta(2,3)**2  ! check
     &  )

        alpha_sq = (beta(1,1) + beta(2,2) + beta(3,3)) / delta_sq

c---------------------------
        ! Finally compute viscosity
        if (alpha_sq .eq. 0.) then
            ediff(i,1,1,e) = param(2)
        else
            ediff(i,1,1,e) = (
     &        param(2) + (
     &          C_vreman * sqrt(max(snrm(i, e), 0.) / alpha_sq)
     &        )
     &      )

        end if
        vdiff(i,1,1,e,1) = ediff(i,1,1,e)
      end do

      return
      end
