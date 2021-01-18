c234567

define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(INTEGER,`integer')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Integrate paths that do not involve level set functions using an c
c     explicit backward in time midpoint method.                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate_paths_midpoint(path, path_gcw, un_0, un_1,
     &         un_2, un_gcw, uh_0, uh_1, uh_2, uh_gcw, dt, dx,
     &         ilower0, ilower1, ilower2, iupper0, iupper1, iupper2)
        implicit none

        INTEGER ilower0, ilower1, ilower2
        INTEGER iupper0, iupper1, iupper2

        INTEGER path_gcw
        REAL path(CELL3d(ilower,iupper,path_gcw),0:2)

        INTEGER un_gcw
        REAL path(SIDE3d0(ilower,iupper,un_gcw))
        REAL path(SIDE3d1(ilower,iupper,un_gcw))
        REAL path(SIDE3d2(ilower,iupper,un_gcw))

        INTEGER uh_gcw
        REAL uh_0(SIDE3d0(ilower,iupper,uh_gcw))
        REAL uh_1(SIDE3d1(ilower,iupper,uh_gcw))
        REAL uh_2(SIDE3d2(ilower,iupper,uh_gcw))

        REAL dt, dx(0:2)

        INTEGER i0, i1, i2
        REAL ux, uy, uz
        REAL xcom, ycom, zcom
        REAL xcom_o, ycom_o, zcom_o

        do i0 = ilower0,iupper0
          do i1 = ilower1,iupper1
            do i2 = ilower2,iupper2
              xcom_o = DBLE(i0) + 0.5d0
              ycom_o = DBLE(i1) + 0.5d0
              zcom_o = DBLE(i2) + 0.5d0
              call find_velocity(i0, i1, i2, un_0, un_1, un_2, un_gcw,
     &            ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &            xcom_o, ycom_o, zcom_o, ux, uy, uz)
              xcom = xcom_o - 0.5d0 * dt * ux / dx(0)
              ycom = ycom_o - 0.5d0 * dt * uy / dx(1)
              zcom = zcom_o - 0.5d0 * dt * uz / dx(2)

              call find_velocity(i0, i1, i2, uh_0, uh_1, uh_2, uh_gcw,
     &            ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &            xcom, ycom, zcom, ux, uy, uz)
              path(i0, i1, i2, 0) = xcom_o - dt * ux / dx(0)
              path(i0, i1, i2, 1) = ycom_o - dt * uy / dx(1)
              path(i0, i1, i2, 2) = zcom_o - dt * uz / dx(2)
            enddo
          enddo
        enddo
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Integrate paths that do not involve level set functions using an c
c     explicit backward in time Euler method.                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate_paths_forward(path, path_gcw, u_0, u_1, u_2,
     &   u_gcw, dt, dx, ilower0, ilower1, ilower2,
     &   iupper0, iupper1, iupper2)
        implicit none

        INTEGER ilower0, ilower1, ilower2
        INTEGER iupper0, iupper1, iupper2

        INTEGER path_gcw
        REAL path(CELL3d(ilower,iupper,path_gcw),0:2)

        INTEGER u_gcw
        REAL u_0(SIDE3d0(ilower,iupper,u_gcw))
        REAL u_1(SIDE3d1(ilower,iupper,u_gcw))
        REAL u_2(SIDE3d2(ilower,iupper,u_gcw))

        REAL dt, dx(0:2)
        INTEGER i0, i1, i2
        REAL ux, uy, uz
        REAL xcom, ycom, zcom

        do i0 = ilower0,iupper0
          do i1 = ilower1,iupper1
            do i2 = ilower2,iupper2
              xcom = DBLE(i0) + 0.5d0
              ycom = DBLE(i1) + 0.5d0
              zcom = DBLE(i2) + 0.5d0
              call find_velocity(i0, i1, i2, u_0, u_1, u_2, u_gcw,
     &            ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &            xcom, ycom, zcom, ux, uy, uz)
              path(i0, i1, i2, 0) = xcom - dt * ux / dx(0)
              path(i0, i1, i2, 1) = ycom - dt * uy / dx(1)
              path(i0, i1, i2, 2) = zcom - dt * uz / dx(2)
            enddo
          enddo
        enddo
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Integrate paths that utilize a level set functions using an      c
c     explicit backward in time midpoint method.                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate_paths_ls_midpoint(path, path_gcw, un_0, un_1,
     &         un_2, un_gcw, uh_0, uh_1, uh_2, uh_gcw, c_data, c_gcw,
     &         dt, dx, ilower0, ilower1, ilower2, iupper0, iupper1,
     &         iupper2)
        implicit none
       
        INTEGER ilower0, ilower1, ilower2
        INTEGER iupper0, iupper1, iupper2
         
        INTEGER path_gcw
        REAL path(CELL3d(ilower,iupper,path_gcw),0:2)
     
        INTEGER un_gcw
        REAL un_0(SIDE3d0(ilower,iupper,un_gcw))
        REAL un_1(SIDE3d1(ilower,iupper,un_gcw))
        REAL un_2(SIDE3d2(ilower,iupper,un_gcw))
     
        INTEGER uh_gcw
        REAL uh_0(SIDE3d0(ilower,iupper,uh_gcw))
        REAL uh_1(SIDE3d1(ilower,iupper,uh_gcw))
        REAL uh_2(SIDE3d2(ilower,iupper,uh_gcw))

        INTEGER c_gcw
        REAL c_data(CELL3d(ilower,iupper,c_gcw),0:2)

        REAL dt, dx(0:2)
        INTEGER i0, i1, i2
        REAL ux, uy, uz
        REAL xcom, ycom, zcom
        REAL xcom_o, ycom_o, zcom_o
         
        do i0 = ilower0,iupper0
          do i1 = ilower1,iupper1
            do i2 = ilower2,iupper2
              xcom_o = c_data(i0, i1, i2, 0)
              ycom_o = c_data(i0, i1, i2, 1)
              zcom_o = c_data(i0, i1, i2, 2)
              call find_velocity(i0, i1, i2, un_0, un_1, un_2, un_gcw,
     &            ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &            xcom_o, ycom_o, zcom_o, ux, uy, uz)
              xcom = xcom_o - 0.5d0 * dt * ux / dx(0)
              ycom = ycom_o - 0.5d0 * dt * uy / dx(1)
              zcom = zcom_o - 0.5d0 * dt * uz / dx(2)

              call find_velocity(i0, i1, i2, uh_0, uh_1, uh_2, uh_gcw,
     &            ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &            xcom, ycom, zcom, ux, uy, uz)
              path(i0, i1, i2, 0) = xcom_o - dt * ux / dx(0)
              path(i0, i1, i2, 1) = ycom_o - dt * uy / dx(1)
              path(i0, i1, i2, 2) = zcom_o - dt * uz / dx(2)
            enddo
          enddo
        enddo
      end
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Integrate paths that utilize a level set functions using an      c
c     explicit backward in time Euler method.                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate_paths_ls_forward(path, path_gcw, u_0, u_1,
     &   u_2, u_gcw, c_data, c_gcw, dt, dx,
     &   ilower0, ilower1, ilower2, iupper0, iupper1, iupper2)
        implicit none

        INTEGER ilower0, ilower1, ilower2
        INTEGER iupper0, iupper1, iupper2
         
        INTEGER path_gcw
        REAL path(CELL3d(ilower,iupper,path_gcw),0:2)
        
        INTEGER u_gcw
        REAL u_0(SIDE3d0(ilower,iupper,u_gcw))
        REAL u_1(SIDE3d1(ilower,iupper,u_gcw))
        REAL u_2(SIDE3d2(ilower,iupper,u_gcw))

        INTEGER c_gcw
        REAL c_data(CELL3d(ilower,iupper,c_gcw),0:2)

        REAL dt, dx(0:2)
        INTEGER i0, i1, i2
        REAL ux, uy, uz
        REAL xcom(0:2)

        do i0 = ilower0,iupper0
          do i1 = ilower1,iupper1
            do i2 = ilower2,iupper2
              xcom(0) = c_data(i0, i1, i2, 0)
              xcom(1) = c_data(i0, i1, i2, 1)
              xcom(2) = c_data(i0, i1, i2, 2)
              call find_velocity(i0, i1, i2, u_0, u_1, u_2, u_gcw,
     &            ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &            xcom(0), xcom(1), xcom(2), ux, uy, uz)
              path(i0, i1, i2, 0) = xcom(0) - dt * ux / dx(0)
              path(i0, i1, i2, 1) = xcom(1) - dt * uy / dx(1)
              path(i0, i1, i2, 2) = xcom(2) - dt * uz / dx(2)
            enddo
          enddo
        enddo
      end
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Interpolate a velocity field to a point (x0, x1) using a         c
c     bilinear interpolant.                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_velocity(i0, i1, i2, u0, u1, u2, u_gcw,
     &      ilower0, ilower1, ilower2, iupper0, iupper1, iupper2,
     &      x0, x1, x2, u0_ret, u1_ret, u2_ret)
        implicit none
        INTEGER i0, i1, i2
        INTEGER ilower0, ilower1, ilower2
        INTEGER iupper0, iupper1, iupper2

        INTEGER u_gcw
        REAL u0(SIDE3d0(ilower,iupper,u_gcw))
        REAL u1(SIDE3d1(ilower,iupper,u_gcw))
        REAL u2(SIDE3d2(ilower,iupper,u_gcw))

        REAL x0, x1, x2

        REAL u0_ret, u1_ret, u2_ret

        REAL xlow, ylow, zlow
        REAL coefs(1:8)
        INTEGER i
        INTEGER j
        INTEGER k
     
c       X coordinate
        if(x1 .gt. (DBLE(i1) + 0.5d0)) then
          if (x2 .gt.(DBLE(i2) + 0.5d0)) then
            i = i0
            j = i1
            k = i2
            xlow = DBLE(i0)
            ylow = DBLE(i1) + 0.5d0
            zlow = DBLE(i2) + 0.5d0
          else
            i = i0
            j = i1
            k = i2-1
            xlow = DBLE(i0)
            ylow = DBLE(i1) + 0.5d0
            zlow = DBLE(i2) - 0.5d0
          endif
        else
          if (x2 .gt. (DBLE(i2) + 0.5d0)) then
            i = i0
            j = i1 - 1
            k = i2
            xlow = DBLE(i0)
            ylow = DBLE(i1) - 0.5d0
            zlow = DBLE(i2) + 0.5d0
          else
            i = i0
            j = i1 - 1
            k = i2 - 1
            xlow = DBLE(i0)
            ylow = DBLE(i1) - 0.5d0
            zlow = DBLE(i2) - 0.5d0
          endif
        endif
        coefs(1) = u0(i,j,k)
        coefs(2) = u0(i+1,j,k) - coefs(1)
        coefs(3) = u0(i,j+1,k) - coefs(1)
        coefs(4) = u0(i,j,k+1) - coefs(1)
        coefs(5) = u0(i+1,j+1,k) - coefs(1) - coefs(2)
        coefs(6) = u0(i,j+1,k+1) - coefs(1) - coefs(3)
        coefs(7) = u0(i+1,j,k+1) - coefs(1) - coefs(4)
        coefs(8) = u0(i+1,j+1,k+1) - SUM(coefs(1:7))
        u0_ret = coefs(1) + coefs(2) * (x0 - xlow)
     &       + coefs(3) * (x1 - ylow) + coefs(4) * (x2 - zlow)
     &       + coefs(5) * (x0 - xlow) * (x1 - ylow)
     &       + coefs(6) * (x1 - ylow) * (x2 - zlow)
     &       + coefs(7) * (x0 - xlow) * (x2 - zlow)
     &       + coefs(8) * (x0 - xlow) * (x1 - ylow) * (x2 - zlow)
c       Y coordinate
        if(x0 .gt. (DBLE(i1) + 0.5d0)) then
          if (x2 .gt.(DBLE(i2) + 0.5d0)) then
            i = i0
            j = i1
            k = i2
            xlow = DBLE(i0) + 0.5d0
            ylow = DBLE(i1)
            zlow = DBLE(i2) + 0.5d0
          else
            i = i0
            j = i1
            k = i2-1
            xlow = DBLE(i0) + 0.5d0
            ylow = DBLE(i1)
            zlow = DBLE(i2) - 0.5d0
          endif
        else
          if (x2 .gt. (DBLE(i2) + 0.5d0)) then
            i = i0 - 1
            j = i1
            k = i2
            xlow = DBLE(i0) - 0.5d0
            ylow = DBLE(i1)
            zlow = DBLE(i2) + 0.5d0
          else
            i = i0 - 1
            j = i1
            k = i2 - 1
            xlow = DBLE(i0) - 0.5d0
            ylow = DBLE(i1)
            zlow = DBLE(i2) - 0.5d0
          endif
        endif
        coefs(1) = u1(i,j,k)
        coefs(2) = u1(i+1,j,k) - coefs(1)
        coefs(3) = u1(i,j+1,k) - coefs(1)
        coefs(4) = u1(i,j,k+1) - coefs(1)
        coefs(5) = u1(i+1,j+1,k) - coefs(1) - coefs(2)
        coefs(6) = u1(i,j+1,k+1) - coefs(1) - coefs(3)
        coefs(7) = u1(i+1,j,k+1) - coefs(1) - coefs(4)
        coefs(8) = u1(i+1,j+1,k+1) - SUM(coefs(1:7))
        u1_ret = coefs(1) + coefs(2) * (x0 - xlow)
     &       + coefs(3) * (x1 - ylow) + coefs(4) * (x2 - zlow)
     &       + coefs(5) * (x0 - xlow) * (x1 - ylow)
     &       + coefs(6) * (x1 - ylow) * (x2 - zlow)
     &       + coefs(7) * (x0 - xlow) * (x2 - zlow)
     &       + coefs(8) * (x0 - xlow) * (x1 - ylow) * (x2 - zlow)
c       Z coordinate
        if(x0 .gt. (DBLE(i1) + 0.5d0)) then
          if (x1 .gt.(DBLE(i2) + 0.5d0)) then
            i = i0
            j = i1
            k = i2
            xlow = DBLE(i0) + 0.5d0
            ylow = DBLE(i1) + 0.5d0
            zlow = DBLE(i2)
          else
            i = i0
            j = i1-1
            k = i2
            xlow = DBLE(i0) + 0.5d0
            ylow = DBLE(i1) - 0.5d0
            zlow = DBLE(i2)
          endif
        else
          if (x1 .gt. (DBLE(i2) + 0.5d0)) then
            i = i0 - 1
            j = i1
            k = i2
            xlow = DBLE(i0) - 0.5d0
            ylow = DBLE(i1) + 0.5d0
            zlow = DBLE(i2)
          else
            i = i0 - 1
            j = i1 - 1
            k = i2
            xlow = DBLE(i0) - 0.5d0
            ylow = DBLE(i1) - 0.5d0
            zlow = DBLE(i2)
          endif
        endif
        coefs(1) = u2(i,j,k)
        coefs(2) = u2(i+1,j,k) - coefs(1)
        coefs(3) = u2(i,j+1,k) - coefs(1)
        coefs(4) = u2(i,j,k+1) - coefs(1)
        coefs(5) = u2(i+1,j+1,k) - coefs(1) - coefs(2)
        coefs(6) = u2(i,j+1,k+1) - coefs(1) - coefs(3)
        coefs(7) = u2(i+1,j,k+1) - coefs(1) - coefs(4)
        coefs(8) = u2(i+1,j+1,k+1) - SUM(coefs(1:7))
        u2_ret = coefs(1) + coefs(2) * (x0 - xlow)
     &       + coefs(3) * (x1 - ylow) + coefs(4) * (x2 - zlow)
     &       + coefs(5) * (x0 - xlow) * (x1 - ylow)
     &       + coefs(6) * (x1 - ylow) * (x2 - zlow)
     &       + coefs(7) * (x0 - xlow) * (x2 - zlow)
     &       + coefs(8) * (x0 - xlow) * (x1 - ylow) * (x2 - zlow)
      end
