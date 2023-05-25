c-----------------------------------------------------------------------
c> Compute \f$\Delta^2\f$, the grid spacing used in a SGS model.
c> @callgraph @callergraph
      subroutine set_grid_spacing
      implicit none

      include 'SIZE'  ! nx1, ny1, nz1, nelev,
      include 'GEOM'  ! xm1, ym1, zm1
      include 'SGS'  ! dg2

      integer e
      integer n, i, j, k, im, ip, jm, jp, km, kp
      real di, dj, dk, ndim_inv

      ndim_inv = 1./ndim

      n = nx1*ny1*nz1*nelv
      call rone(dg2,n)
c     return               ! Comment this line for a non-trivial Delta defn
      do e=1,nelv

         do k=1,nz1
           km = max(1  ,k-1)
           kp = min(nz1,k+1)

           do j=1,ny1
             jm = max(1  ,j-1)
             jp = min(ny1,j+1)

             do i=1,nx1
               im = max(1  ,i-1)
               ip = min(nx1,i+1)

               di = (xm1(ip,j,k,e)-xm1(im,j,k,e))**2
     $            + (ym1(ip,j,k,e)-ym1(im,j,k,e))**2
     $            + (zm1(ip,j,k,e)-zm1(im,j,k,e))**2

               dj = (xm1(i,jp,k,e)-xm1(i,jm,k,e))**2
     $            + (ym1(i,jp,k,e)-ym1(i,jm,k,e))**2
     $            + (zm1(i,jp,k,e)-zm1(i,jm,k,e))**2

               dk = (xm1(i,j,kp,e)-xm1(i,j,km,e))**2
     $            + (ym1(i,j,kp,e)-ym1(i,j,km,e))**2
     $            + (zm1(i,j,kp,e)-zm1(i,j,km,e))**2
c               write(6,*) ip,im,jp,jm
               di = sqrt(di)/(ip-im)
               dj = sqrt(dj)/(jp-jm)
               dk = sqrt(dk)/(kp-km)
               dg2(i,j,k,e) = (di*dj*dk)**(2*ndim_inv)

             enddo
           enddo
         enddo
      enddo

      call dsavg(dg2)  ! average neighboring elements

      do e=1,nelv
        dg2_max(e) = maxval(dg2(:,:,:,e))
      enddo

#ifdef DEBUG
      if (nid.eq.0) print *, "set_grid_spacing :dg2_max =", dg2_max
      print *, "set_grid_spacing: rank =", nid,
     &     ",local max =", maxval(dg2),
     &     ",local min =", minval(dg2)
      if (nid.eq.0) then
        print *, "set_grid_spacing; lx1, ly1, lz1, lelv: ",
     &      lx1,ly1,lz1,lelv
        open(unit=5, file="set_grid_spacing.dat", form="unformatted")
        write(5) xm1
        write(5) ym1
        write(5) zm1
        write(5) dg2
        write(5) dg2_max
        close(5)
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
c> Setup test filter
c> @callgraph @callergraph
      subroutine set_ds_filt(fh,fht,nt,diag,nx)
      implicit none

      INCLUDE 'SIZE'

      integer nt, nx
      real fh(nx*nx),fht(nx*nx),diag(nx)

c Construct transfer function
      call rone(diag,nx)

c      diag(nx-0) = 0.01
c      diag(nx-1) = 0.10
c      diag(nx-2) = 0.50
c      diag(nx-3) = 0.90
c      diag(nx-4) = 0.99
c      nt = nx - 2

      diag(nx-0) = 0.05
      diag(nx-1) = 0.50
      diag(nx-2) = 0.95
      nt = nx - 1

      call build_1d_filt(fh,fht,diag,nx,nid)

      return
      end
c-----------------------------------------------------------------------
