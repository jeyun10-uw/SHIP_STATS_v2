
subroutine statistics()

use vars
use grid, only: doShipTrackConditionals, day, nstat, nstep, nstatfrq
use micro_params, only: shipv2_time0, NC_env, NAc_thresh
use rad, only: qrad, do_output_clearsky_heating_profiles, radqrclw, radqrcsw, &
               swUp3D, swDown3D, lwUp3D, lwDown3D
use microphysics, only: Get_dryaerosol, Get_qt, Get_cloudaerosol, Get_accumaerosol 
use domain, only: nx_gl
use tracers
use params
use hbuffer
use instrument_diagnostics, only: compute_instr_diags
implicit none	
	
	real mse(nzm)
	real dse(nzm)
	real sse(nzm)
	real tpz(nzm)
	real tlz(nzm)
	real tvz(nzm)
        real qcz(nzm)
	real qiz(nzm)
	real tez(nzm)
	real qvz(nzm)
	real qrz(nzm)
	real qsz(nzm)
	real relhz(nzm)

	real u2z(nzm)
	real v2z(nzm)
	real w2z(nzm)
	real w22(nzm)
	real w3z(nzm)
	real skw(nzm)
	real t2z(nzm)
	real tqz(nzm)
	real q2z(nzm)
	real qc2z(nzm)
	real qi2z(nzm)
	real qs2z(nzm)
	real tkez(nzm)
	real fadv(nz)
	real shear(nz)
	real shearx(nzm)
	real sheary(nzm)
	real presx(nzm)
	real presy(nzm)
	real twgrad(nzm)
	real qwgrad(nzm)
	real swgrad(nzm)
	
	
	real tvwle(nzm)
	real qcwle(nzm)
	real qiwle(nzm)
	real aup(nzm)
	real wcl(nzm)
	real ucl(nzm)
	real vcl(nzm)
	real tcl(nzm)
	real tacl(nzm)
	real tvcl(nzm)
	real qcll(nzm)
	real qccl(nzm)
	real qicl(nzm)
	real qpcl(nzm)
	real twcl(nzm)
	real qwcl(nzm)
	real tvwcl(nzm)
	real qcwcl(nzm)
	real qiwcl(nzm)
	real wacl(nzm)	
	real cld(nzm)
 	real cldd(nzm)
	real cldcumup1(nzm), cldcumup2(nzm), tmpcwp
        real cldcumdn1(nzm), cldcumdn2(nzm), cwp_threshold1, cwp_threshold2
	real hydro(nzm)
	real qsatwz(nzm)

	real tvirt(nx,ny,nzm)
	
	integer i,j,k,n,ntr
	real qcc,qii,qrr,qss,lstarn,lstarp,coef,coef1
	real factor_xy, factor_n, tmp(4), tmp1(4)
        real buffer(nzm,6),buffer1(nzm,6)
	real prof1(nzm),prof2(nzm),prof3(nzm),prof4(nzm)	
	real cwpmax,cwp(nx,ny),cwpl(nx,ny),cwpm(nx,ny),cwph(nx,ny)
	logical condition, condition_cl, condition_cloud, condition_track
	real zero(nzm)

	integer topind(nx,ny),z_inv_ind(nx,ny),z_base_ind(nx,ny),z_top_ind(nx,ny),ncloud	
        real zzz,grad_max(nx,ny),grad,tmpqcl,tmpqci, qclz(nzm), qciz(nzm)

!========================================================================
! UW ADDITIONS

 real tvcla(nzm)!kzm added Apr. 7,2004 for thetav anomalies
 real wstar3(nzm) !bloss added 11/04/05 

 !----------------------SHIP TRACK CONDITIONALS (mcmichael) -----------------------
 !radiation streams sampled conditionally
 real, dimension(nzm) :: sw_u, sw_d, lw_u, lw_d

 !variances with ship track conditionals
 !qs2 -> total water variance
 !ts2 -> l.s.e variance
 !ws2 -> vertical velocity variance
 !us2, vs2 -> horizontal velocity variances
 real, dimension(nzm) :: qs2, ts2, ws2 
 real, dimension(nzm) :: us2, vs2

 !vertical turbulent fluxes
 !wtl -> l.s.e vertical flux
 !wqt -> vertical total water flux
 !wtv -> vertical buoyancy flux
 real, dimension(nzm) :: wtl, wqt, wtv

 !liquid precipitation flux for ship tracks
 real, dimension(nzm) :: pfs
 
 !Ship track identification variables
 real :: std_t = 3.0 !standard deviation 
 !moved aero_col, std_aero, aero_thresh, stats_flag to vars.f90
 real :: buffer30(5), buffer31(5), NC_env2, NAC_mean2, NAC2_mean, NAC_std2, &
         ncenv_numerator, ncenv_denominator, nacenv_denominator, nacenv_numerator, nacenv2_numerator, &
         rhodz
 real, dimension(nx) :: edge_factor
 real, dimension(nx,ny) :: aero_xy, aero_xy_col

 !column-by-column MLM stuff
 real, dimension(nx,ny,nzm) :: T_shift !correcting tabs grid 
 real, dimension(nx,ny,nzm) :: T_diff !tabs grid for derivative
 real, dimension(nx,ny) :: height_inv !array of inversion base heights
 real, dimension(nx,ny) :: height_cb, height_ct, height_cb1, height_ct1, lwp, cl_depth, cb_h, NCenv, NAcenv
 real, dimension(nx,ny) :: tlcl, qtcl, tlsc, qtsc, uclw, vclw, uscw, vscw !cl, subcloud layer averages
 real, dimension(nx,ny) :: qt_w, tl_w, u_w, v_w !weighted BL column means
 real, dimension(nzm) :: z_grid, z_grid1, z_diff1 !correcting z grid
 real, dimension(nzm) :: zb, tcb, zbct !inversion base height
 real, dimension(nzm) :: wlo, wloct !local vertical velocity
 real, dimension(nzm) :: wla, wlact !large-scale subsidence
 real, dimension(nzm) :: swu0, swd0, lwu0, lwd0, swut, swdt, lwut, lwdt, &
                         swutct, swdtct, lwutct, lwdtct !rad terms
 real, dimension(nzm) :: rdiv, ldiv, sdiv, rdivct, sdivct, ldivct !radiative divergence terms
 real, dimension(nzm) :: sh_f, lh_f !surface flux terms 
 real, dimension(nzm) :: pr_s, pr_z, pr_zct, pr_zb !surface, inversion top, cb prec. flux
 real, dimension(nzm) :: cb_h_avg, cl_d_avg, lwp_avg, &
                         tlcl_avg, tlsc_avg, qtcl_avg, qtsc_avg, &
                         tl_avg, qt_avg
 
 !horizontal advection related variables
 real, dimension(nzm) :: u_adv_t, u_adv_q, v_adv_t, v_adv_q, w_adv_t, w_adv_q
 real, parameter :: der_thresh = 0.03 !derivative threshold for dT/dz
 integer, parameter :: h_cb = 10 !starting index for k loop (~cloud-base)
 real :: der_temp, z_tot, rho_tot, rho_tot2, rho_tot3, lwc_tot, z_fraction, rho_fraction, weight
 
 !relaxation timescale variables
 real, dimension(nzm) :: difftke, disstke, restke

 !dry aerosol concentration (#/mg) passed from microphysics.f90
 real, dimension(nx,ny,nzm) :: Dryaero

 !total water passed from microphysics.f90 with buffer points
 real, dimension(0:nx+1,1-YES3D:ny+YES3D,nzm) :: tot_w

 !cloud number aerosol concentration (#/mg) passed from microphysics.f90
 real, dimension(nx,ny,nzm) :: Cloudaero

 !accumulation mode aerosol concentration (#/mg) passed from microphysics.f90
 real, dimension(nx,ny,nzm) :: Accumaero

 !variables associated with environmental aerosol concentration calculations
 real :: NAc_env, std_nac_env
 real, dimension(nx) :: tmpx_pos

 !fraction of grid on edges of the domain used to calculate environmental NC
 real, parameter :: edge_frac = 0.125 

 !-------------------------------------------------------------------------------

 !bloss: momentum flux statistics for cloud, up/downdraft cores
 real, dimension(nzm) :: uwsbcl, vwsbcl, uwlecl, vwlecl, &
      uadv, vadv, udiff, vdiff
 real :: uwsubgrid, vwsubgrid, uwresolved, vwresolved

 !bloss: new stuff for conditionally-averaged statistics (i.e. cloud, core, etc.)
 integer ncond, jb, kc, kb, it, jt 
 character(LEN=6) :: statname
 real :: tmprhow, tmpmse, tmpqt

 !bloss: conditional u,v anomalies, pressure gradients
 real, dimension(nzm) :: ucla, vcla, dpdxcl, dpdycl, dpdzcl
               
 !bloss: frozen moist static energy
 real, dimension(nzm) :: fmse, fmsecla

 !bloss: mass flux and mass-flux weighted stats in conditional category
real, dimension(nzm) :: rhowcl, rhowmsecl, rhowtlcl, rhowqtcl,  &
     rhowmsecla, rhowtlcla, rhowqtcla, rhouwcl, rhovwcl, rhowwcl, &
     rhowtvcl, rhowtvcla

real :: relhobs(nzm)

! END UW ADDITIONS
!========================================================================

        call t_startf('statistics')

	factor_xy = 1./float(nx*ny)
	factor_n = 1./float(nsubdomains)
	
!bloss: Additional calls to boundaries so that clean momentum flux
! budgets can be computed. 
!----------------------------------------------------------
!      Update the subdomain's boundaries for velocity

     call boundaries(0)

!---------------------------------------------------------
!        Update boundaries for the SGS exchange coefficients:

     call boundaries(4)

!---------------------------------------------------------
!mcmichael: Update boundaries for scalar quantities for 
!advection calculations

     call boundaries(2)

     !grab total water from microphysics.f90
     tot_w(0:nx+1,1-YES3D:ny+YES3D,:) = Get_qt()
!-----------------------------------------------
!	Mean thermodynamics profiles:
!-----------------------------------------------	
		
	do k=1,nzm
	 dse(k)=0.
	 mse(k)=0.
	 sse(k)=0.
	 tpz(k) = 0.
	 tlz(k) = 0.
	 tvz(k) = 0.
	 tez(k) = 0.
	 qvz(k) = 0.
	 qcz(k) = 0.
	 qiz(k) = 0.
	 qrz(k) = 0.
	 qsz(k) = 0.
	 qsatwz(k)=0.
	 relhz(k)=0.
	 prof1(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
	 zero(k)=0.
	 do j=1,ny
	  do i=1,nx
	   qcc=qcl(i,j,k)
	   qii=qci(i,j,k)
	   qrr=qpl(i,j,k)
	   qss=qpi(i,j,k)
           qrz(k)=qrz(k)+qrr
           qsz(k)=qsz(k)+qss
           qcz(k)=qcz(k)+qcc
           qiz(k)=qiz(k)+qii
	   prof1(k)=prof1(k)+qcc+qii
	   prof2(k)=prof2(k)+qrr+qss
	   prof3(k)=prof3(k)+qcc+qii+qrr+qss
	   tmp(1)=tabs(i,j,k)*prespot(k)
           tpz(k)=tpz(k)+tmp(1)
	   tlz(k)=tlz(k)+tmp(1)*(1.-fac_cond*(qcl(i,j,k)+qci(i,j,k))/tabs(i,j,k))
           tvirt(i,j,k)=tmp(1)*(1.+epsv*qv(i,j,k)-(qcl(i,j,k)+qci(i,j,k))-(qpl(i,j,k)+qpi(i,j,k)))
	   tvz(k)=tvz(k)+tvirt(i,j,k)
           tez(k)=tez(k)+tabs(i,j,k)+gamaz(k)+fac_cond*qv(i,j,k)-fac_fus*(qii+qss)
	   qvz(k) =qvz(k)+qv(i,j,k)
     	   dse(k)=dse(k)+tabs(i,j,k)+gamaz(k)
	   mse(k)=mse(k)+tabs(i,j,k)+gamaz(k)+fac_cond*qv(i,j,k)
	   sse(k)=sse(k)+tabs(i,j,k)+gamaz(k)+fac_cond*qsatw(tabs(i,j,k),pres(k))
	   qsatwz(k) = qsatwz(k)+qsatw(tabs(i,j,k),pres(k))
	   relhz(k)=relhz(k)+qv(i,j,k)/qsatw(tabs(i,j,k),pres(k))
	  end do
	 end do
	end do	
	

	call hbuf_avg_put('TL',t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
	call hbuf_avg_put('TABS',tabs,1,nx, 1,ny, nzm,1.)
	call hbuf_avg_put('U',u+ug,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1.)
	call hbuf_avg_put('V',v+vg,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1.)
	call hbuf_avg_put('QT',qv+qcl+qci,1,nx,1,ny,nzm,1.e3)
	call hbuf_put('TABSOBS',tg0,1.)
	call hbuf_put('QVOBS',qg0,1.e3)
        call hbuf_avg_put('UXGRID',u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,1.)
	call hbuf_avg_put('VXGRID',v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,1.)

        if(nudge_to_sounding_winds) then
          !bloss: If nudging to sounding wind, output sounding as [UV]OBS
          call hbuf_put('UOBS',usounding0+ug,1.)
          call hbuf_put('VOBS',vsounding0+vg,1.)
          if(donudging_uv) then
            call hbuf_put('UBIAS',u0-usounding0,1.)
            call hbuf_put('VBIAS',v0-vsounding0,1.)
          end if
        else
          !bloss: Otherwise, output geostrophic wind as [UV]OBS
          call hbuf_put('UOBS',ug0+ug,1.)
          call hbuf_put('VOBS',vg0+vg,1.)
          if(donudging_uv) then
            call hbuf_put('UBIAS',u0-ug0,1.)
            call hbuf_put('VBIAS',v0-vg0,1.)
          end if
        end if

        if(docoriolis) then
          !bloss: if coriolis is on, output geostrophic wind explicitly
          call hbuf_put('UGEOSTR',ug0+ug,1.)
          call hbuf_put('VGEOSTR',vg0+vg,1.)
        end if

        call hbuf_put('WOBS',wsub,1.)
	call hbuf_put('TTEND',ttend,86400.)
	call hbuf_put('QTEND',qtend,86400.*1.e3)

	call hbuf_put('DSE',dse,factor_xy)
	call hbuf_put('MSE',mse,factor_xy)
	call hbuf_put('SSE',sse,factor_xy)
	call hbuf_put('THETA',tpz,factor_xy)
	call hbuf_put('THETAL',tlz,factor_xy)
	call hbuf_put('THETAV',tvz,factor_xy)
	call hbuf_put('THETAE',tez,factor_xy)

	call hbuf_put('PRES',pres,1.)
	call hbuf_put('RHO',rho,1.)
	call hbuf_put('QV',qvz,1.e3*factor_xy)
	call hbuf_put('QCL',qcz,1.e3*factor_xy)
	call hbuf_put('QCI',qiz,1.e3*factor_xy)
	call hbuf_put('QPL',qrz,1.e3*factor_xy)
	call hbuf_put('QPI',qsz,1.e3*factor_xy)
	call hbuf_put('QN',prof1,1.e3*factor_xy)
	call hbuf_put('QP',prof2,1.e3*factor_xy)
	call hbuf_put('QCOND',prof3,1.e3*factor_xy)
	call hbuf_put('QSAT',qsatwz,1.e3*factor_xy)
	call hbuf_put('RELH',relhz,100.*factor_xy)

 !bloss(2018-11-29): Add bias outputs
 do k = 1,nzm
   relhobs(k)=qg0(k)/qsatw(tg0(k),pres(k))
 end do
 call hbuf_put('TBIAS',tabs0-tg0,1.)
 call hbuf_put('QBIAS',factor_xy*(qvz+qcz)-qg0,1.e3)
 call hbuf_put('RELHBIAS',factor_xy*relhz-relhobs,100.)

!-------------------------------------------------------------
!	Fluxes:
!-------------------------------------------------------------

	do k=1,nzm
	  tmp(1) = dz/rhow(k)
	  tmp(2) = tmp(1) / dtn
	  uwsb(k) = uwsb(k) * tmp(1)
	  vwsb(k) = vwsb(k) * tmp(1)
	  twsb(k) = twsb(k) * tmp(1) * rhow(k) * cp
	  uwle(k) = uwle(k)*tmp(1) + uwsb(k)
	  vwle(k) = vwle(k)*tmp(1) + vwsb(k)
	  twle(k) = twle(k)*tmp(2)*rhow(k)*cp + twsb(k)
	  if(dotracers) then
           do ntr=1,ntracers
	    trwsb(k,ntr) = trwsb(k,ntr) * tmp(1)*rhow(k)
     	    trwle(k,ntr) = trwle(k,ntr) * tmp(2)*rhow(k) + trwsb(k,ntr)
           end do
          end if
	end do
	uwle(nz) = 0.
	vwle(nz) = 0.
	uwsb(nz) = 0.
	vwsb(nz) = 0.

	call hbuf_put('UW',uwle,factor_xy)
	call hbuf_put('VW',vwle,factor_xy)
	call hbuf_put('UWSB',uwsb,factor_xy)
	call hbuf_put('VWSB',vwsb,factor_xy)
	call hbuf_put('TLFLUX',twle,factor_xy)
	call hbuf_put('TLFLUXS',twsb,factor_xy)
	call hbuf_put('PRECIP',precflux,factor_xy/dt*dz*86400./(nstatis+1.e-5))
	
	do j=1,ny
	 do i=1,nx
	  precsfc(i,j)=precsfc(i,j)*dz/dt*86400./(nstatis+1.e-5)
          if(precsfc(i,j).gt.0.1) s_ar = s_ar + 1.
          if(precsfc(i,j)/86400./rhow(1).gt.3.65e-5) s_arthr = s_arthr + 1.
	 end do
	end do
        precmax = maxval(precsfc(:,:))
        precmean = precmean+sum(precsfc(:,:))
        prec2 = prec2+sum(precsfc(:,:)**2)


	do k=1,nzm
	 tvz(k) = 0.
	 qcz(k) = 0.
	 qiz(k) = 0.
	 qsatwz(k) = 0.
	 prof1(k)=0.
	 prof2(k)=0.
	 do j=1,ny
	  do i=1,nx
	    tvz(k) = tvz(k) + tvirt(i,j,k)
	    qcz(k) = qcz(k) + qcl(i,j,k)	    
            qiz(k) = qiz(k) + qci(i,j,k)
	    qsatwz(k) = qsatwz(k)+qsatw(tabs(i,j,k),pres(k))
            qrz(k) = qrz(k) + qpl(i,j,k)
	  end do
	 end do
	 tvz(k) = tvz(k)*factor_xy
	 qcz(k) = qcz(k)*factor_xy
	 qiz(k) = qiz(k)*factor_xy	 
	 qsatwz(k) = qsatwz(k)*factor_xy
         qrz(k) = qrz(k)*factor_xy
	end do
	if(dompi) then
	  coef1 = 1./float(nsubdomains)
	  do k=1,nzm
	    buffer(k,1) = tvz(k)
	    buffer(k,2) = qcz(k)
	    buffer(k,3) = qiz(k)
	    buffer(k,4) = qsatwz(k)
            buffer(k,5) = qrz(k)
	  end do
	  call task_sum_real(buffer,buffer1,nzm*5)
	  do k=1,nzm
	    tvz(k) = buffer1(k,1) * coef1
	    qcz(k) = buffer1(k,2) * coef1
	    qiz(k) = buffer1(k,3) * coef1
	    qsatwz(k) = buffer1(k,4) * coef1
            qrz(k) = buffer1(k,5) * coef1
	  end do
	end if ! dompi

	tvwle(1) = 0.
        wstar3(1) = 0. !bloss
	qcwle(1) = 0.
	qiwle(1) = 0.
	do k=2,nzm
	 tvwle(k) = 0.
         wstar3(k) = 0. !bloss
	 qcwle(k) = 0.
	 qiwle(k) = 0.
	 do j=1,ny
	  do i=1,nx
	    tvwle(k) = tvwle(k) + 0.5*w(i,j,k)* &
		(tvirt(i,j,k-1)-tvz(k-1)+tvirt(i,j,k)-tvz(k))
	    qcwle(k) = qcwle(k) + 0.5*w(i,j,k)* &
		(qcl(i,j,k-1)-qcz(k-1)+ qcl(i,j,k)-qcz(k))	  
	    qiwle(k) = qiwle(k) + 0.5*w(i,j,k)* &
                (qci(i,j,k-1)-qiz(k-1)+qci(i,j,k)-qiz(k))
	    prof1(k)=prof1(k)+rho(k)*0.5* &
                (w(i,j,k)**2+w(i,j,k+1)**2)*(t(i,j,k)-t0(k))
          end do
	 end do
         wstar3(k) = wstar3(k-1) + 2.5*dz*adzw(k)*bet(k)*tvwle(k) !bloss
	 tvwle(k) = tvwle(k)*rhow(k)*cp
	 qcwle(k) = qcwle(k)*rhow(k)*lcond
	 qiwle(k) = qiwle(k)*rhow(k)*lcond
	end do	

	call hbuf_put('TVFLUX',tvwle,factor_xy)
	call hbuf_put('QCFLUX',qcwle,factor_xy)
	call hbuf_put('QIFLUX',qiwle,factor_xy)

        !bloss: UW additions
	call hbuf_put('WSTAR3',wstar3,factor_xy) !bloss

!---------------------------------------------------------
!	Mean turbulence related profiles:
!-----------------------------------------------------------


	do k=1,nzm

	 u2z(k) = 0.
	 v2z(k) = 0.
	 w2z(k) = 0.
	 w22(k) = 0.
	 w3z(k) = 0.
	 aup(k) = 0.
	 t2z(k) = 0.
	 tqz(k) = 0.
	 q2z(k) = 0.
	 qc2z(k) = 0.
	 qi2z(k) = 0.
	 qs2z(k) = 0.
	 do j=1,ny
	  do i=1,nx
	    u2z(k) = u2z(k)+(u(i,j,k)-u0(k))**2	  
	    v2z(k) = v2z(k)+(v(i,j,k)-v0(k))**2
	    w2z(k) = w2z(k)+0.5*(w(i,j,k+1)**2+w(i,j,k)**2)
	    w22(k) = w22(k)+w(i,j,k)**2
	    w3z(k) = w3z(k)+0.5*(w(i,j,k+1)**3+w(i,j,k)**3)	  
	    t2z(k) = t2z(k)+(t(i,j,k)-t0(k))**2	  
	    tqz(k) = tqz(k)+(t(i,j,k)-t0(k))*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k))	  
	    q2z(k) = q2z(k)+(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k))**2
	    if(w(i,j,k)+w(i,j,k+1).gt.0) aup(k) = aup(k) + 1
	  end do
	 end do
	 skw(k) = w3z(k)/(w2z(k)*factor_xy+1.e-5)**1.5
	 tkez(k)= 0.5*(u2z(k)+v2z(k)*YES3D+w2z(k))
	 tvwle(k) = tvwle(k) * bet(k) /(rho(k)*cp)
	 do j=1,ny
	  do i=1,nx
	    qc2z(k) = qc2z(k)+(qcl(i,j,k)-qcz(k))**2
	    qi2z(k) = qi2z(k)+(qci(i,j,k)-qiz(k))**2
	    qs2z(k) = qs2z(k)+(qsatw(tabs(i,j,k),pres(k))-qsatwz(k))**2
	  end do
	 end do
	 
	end do

	call hbuf_put('U2',u2z,factor_xy)
	call hbuf_put('V2',v2z,factor_xy)
	call hbuf_put('W2',w2z,factor_xy)
	call hbuf_put('W3',w3z,factor_xy)
	call hbuf_put('WSKEW',skw,factor_xy)
	call hbuf_put('AUP',aup,factor_xy)

	call hbuf_put('TL2',t2z,factor_xy)
	call hbuf_put('TQ',tqz,factor_xy)
	call hbuf_put('QT2',q2z,1.e6*factor_xy)
	call hbuf_put('QC2',qc2z,1.e6*factor_xy)
	call hbuf_put('QI2',qi2z,1.e6*factor_xy)
	call hbuf_put('QS2',qs2z,1.e6*factor_xy)
	
	call hbuf_put('TKE',tkez,factor_xy)
!-----------------------------------------------------------------
!  TKE balance:

        shear(1)=0.
        shear(nz)=0.
        do k=2,nzm
          shear(k)=-( (uwle(k)-uwsb(k))*(u0(k)-u0(k-1)) &
            +(vwle(k)-vwsb(k))*(v0(k)-v0(k-1))*YES3D )*factor_xy /(dz*adzw(k))
        end do
        do k=1,nzm
          shear(k)=0.5*(shear(k)+shear(k+1))
          tkeleadv(k)=tkeleadv(k)-shear(k)
          tkelediff(k)=tkelediff(k)-tkelediss(k)
        end do

        call hbuf_put('ADVTR',tkeleadv,1.)
        call hbuf_put('PRESSTR',tkelepress,1.)
        call hbuf_put('BUOYA',tkelebuoy,1.)
        call hbuf_put('SHEAR',shear,1.)
        call hbuf_put('DISSIP',tkelediss,1.)
        call hbuf_put('DIFTR',tkelediff,1.)
	
        ! damping of TKE by damping layer at top of domain
        call hbuf_put('TKEDAMP',tkedamp,factor_xy) !bloss(2019-04-02)

	fadv(1)=0.
	fadv(nz)=0.
!-----------------------------------------------------------------
!  Momentum flux balance:

!  UW advection d(w'w'u')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                         ( u(i,j,k-1)-u0(k-1)+u(i,j,k)-u0(k))
	  end do
	 end do
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 shearx(k)=momleadv(k,1)-coef
	 momleadv(k,1)=coef
	end do	


!  VW advection d(w'w'v')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &	
                         ( v(i,j,k-1)-v0(k-1)+v(i,j,k)-v0(k))
	  end do
	 end do
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 sheary(k)=momleadv(k,2)-coef
	 momleadv(k,2)=coef
	end do	


!  UW advection d(p'u')/dz:

	do k=1,nz
	 fadv(k)=0.
	 if(k.eq.1) then
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+(1.5*(u(i,j,k)-u0(k))*p(i,j,k)*rho(k)- &
                     0.5*(u(i,j,k+1)-u0(k+1))*p(i,j,k+1)*rho(k+1))
	    end do
	   end do
	 else if(k.eq.nz) then
	   do j=1,ny
	    do i=1,nx
	     fadv(k)=fadv(k)+(1.5*(u(i,j,k-1)-u0(k-1))*p(i,j,k-1)*rho(k-1)- &
                  0.5*(u(i,j,k-2)-u0(k-2))*p(i,j,k-2)*rho(k-2))
	    end do
	   end do
	 else
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+0.5*((u(i,j,k)-u0(k))*p(i,j,k)*rho(k)+ &
                       (u(i,j,k-1)-u0(k-1))*p(i,j,k-1)*rho(k-1))
	    end do
	   end do
	 end if
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 presx(k)=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	end do	


!  VW advection d(p'v')/dz:

	do k=1,nz
	 fadv(k)=0.
	 if(k.eq.1) then
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+(1.5*(v(i,j,k)-v0(k))*p(i,j,k)*rho(k)- & 
                    0.5*(v(i,j,k+1)-v0(k+1))*p(i,j,k+1)*rho(k+1))
	    end do
	   end do
	 else if(k.eq.nz) then
	   do j=1,ny
	    do i=1,nx
             fadv(k)=fadv(k)+(1.5*(v(i,j,k-1)-v0(k-1))*p(i,j,k-1)*rho(k-1)- &
                     0.5*(v(i,j,k-2)-v0(k-2))*p(i,j,k-2)*rho(k-2))
	    end do
	   end do
	 else
	   do j=1,ny
	    do i=1,nx
	      fadv(k)=fadv(k)+0.5*((v(i,j,k)-v0(k))*p(i,j,k)*rho(k)+ &
                       (v(i,j,k-1)-v0(k-1))*p(i,j,k-1)*rho(k-1))
	    end do
	   end do
	 end if
	 fadv(k)=fadv(k)*factor_xy
	end do
	do k=1,nzm
	 presy(k)=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	end do	

	do k=1,nzm
	  momlepress(k,1)=momlepress(k,1)-presx(k)
	  momlepress(k,2)=momlepress(k,2)-presy(k)
	  momlepress(k,3)=momlepress(k,3)-tkelepress(k)
	end do

	call hbuf_put('WUADV',momleadv(1,1),1.)
	call hbuf_put('WUANIZ',momlepress(1,1),1.)
	call hbuf_put('WUBUOY',momlebuoy(1,1),1.)
	call hbuf_put('WUSHEAR',shearx,1.)
	call hbuf_put('WUPRES',presx,1.)
	call hbuf_put('WUDIFF',momlediff(1,1),1.)

	call hbuf_put('WVADV',momleadv(1,2),1.)
	call hbuf_put('WVANIZ',momlepress(1,2),1.)
	call hbuf_put('WVBUOY',momlebuoy(1,2),1.)
	call hbuf_put('WVSHEAR',sheary,1.)
	call hbuf_put('WVPRES',presy,1.)
	call hbuf_put('WVDIFF',momlediff(1,2),1.)

	call hbuf_put('W2BUOY',momlebuoy(1,3),2.)
	call hbuf_put('W2ADV',momleadv(1,3),2.)
	call hbuf_put('W2REDIS',momlepress(1,3),2.)
	call hbuf_put('W2PRES',tkelepress,2.)
	call hbuf_put('W2DIFF',momlediff(1,3),2.)

!-----------------------------------------------------------
! T2 and Q2 variance budget:


	do k=1,nzm
	  q2lediff(k)=q2lediff(k)-q2lediss(k)
	  t2lediff(k)=t2lediff(k)-t2lediss(k)
	end do	
	

	call hbuf_put('T2ADVTR',t2leadv,1.)
	call hbuf_put('T2GRAD',t2legrad,1.)
	call hbuf_put('T2DISSIP',t2lediss,1.)
	call hbuf_put('T2DIFTR',t2lediff,1.)
	call hbuf_put('T2PREC',t2leprec,1.)

	call hbuf_put('Q2ADVTR',q2leadv,1.)
	call hbuf_put('Q2GRAD',q2legrad,1.)
	call hbuf_put('Q2DISSIP',q2lediss,1.)
	call hbuf_put('Q2DIFTR',q2lediff,1.)
	call hbuf_put('Q2PREC',q2leprec,1.)

!------------------------------------------------------------------
! HW and QW budgets:


	fadv(1)=0.
	fadv(nz)=0.


!  HW advection d(w'w'h')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
                 	      ( t(i,j,k-1)-t0(k-1)+t(i,j,k)-t0(k))
	  end do
	 end do
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 twgrad(k)=twleadv(k)-coef
	 twleadv(k)=coef
	end do	


!  QW advection d(w'w'q')/dz:

	do k=2,nzm
	 fadv(k)=0.
	 do j=1,ny
	  do i=1,nx
	    fadv(k)=fadv(k)+w(i,j,k)**2*rhow(k)*0.5* &
       	     ( qv(i,j,k-1)+qcl(i,j,k)+qci(i,j,k)-q0(k-1)+qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k))
	  end do
	 end do
	end do
	do k=1,nzm
	 coef=-(fadv(k+1)-fadv(k))/(adz(k)*dz*rho(k))
	 qwgrad(k)=qwleadv(k)-coef
	 qwleadv(k)=coef
	end do	



	call hbuf_put('TWADV',twleadv,factor_xy)
	call hbuf_put('TWDIFF',twlediff,factor_xy)
	call hbuf_put('TWGRAD',twgrad,factor_xy)
	call hbuf_put('TWBUOY',twlebuoy,factor_xy)
	call hbuf_put('TWPRES',twlepres,factor_xy)
	call hbuf_put('TWPREC',twleprec,factor_xy)

	call hbuf_put('QWADV',qwleadv,factor_xy)
	call hbuf_put('QWDIFF',qwlediff,factor_xy)
	call hbuf_put('QWGRAD',qwgrad,factor_xy)
	call hbuf_put('QWBUOY',qwlebuoy,factor_xy)
	call hbuf_put('QWPRES',qwlepres,factor_xy)
	call hbuf_put('QWPREC',qwleprec,factor_xy)

!-------------------------------------------------------------
!	Conditional statistics:
!-------------------------------------------------------------
!bloss:  Major modification of conditional statistics to allow 
         ! a standard set of averages to be defined using many
         ! different conditions.  Here, the conditions are cloudy,
         ! updraft core, downdraft core, saturated updraft, saturated
         ! downdraft and unsaturated environment.

         ! initialize mask array, which will be one only where the conditional
         !   is satisfied, and zero everywhere else.
         !   this is used in MICRO_M2005 to compute conditional averages of 
         !     microphysical transfer rates
         
!mcmichael: added conditional statistics for ship track stats 
         !conditions are non-cloudy ship track, non-cloudy no ship,
         !cloudy ship track, and cloudy no ship
         !Ship track is determined by first calculating the boundary-layer
         !averaged NAc on the edges of the domain (fraction of edge controlled 
         !by edge_frac = 0.125) at each statistics call. The ship threshold is
         !then calculated as 3 times the standard deviation of the environmental NAc. 
         !--------------------------------------------------------------------------------------
         if(doShipTrackConditionals) then

           aero_xy_col(:,:) = 0.

                !shifting verticl grid to grid box center
                do k = 1,nzm-1
                   z_grid1(k) = (z(k+1) + z(k))/2.0
                   T_shift(:,:,k) = (tabs(:,:,k+1) + tabs(:,:,k))/2.0
                enddo
                
                !compute dz
                do k = 1,nzm
                   if (k.eq.1) then
                      z_diff1(k) = z_grid1(1)
                      T_diff(:,:,k) = 0.
                   else
                      z_diff1(k) = z(k) - z(k-1)
                      T_diff(:,:,k) = T_shift(:,:,k) - T_shift(:,:,k-1)
                   endif
                enddo

                !compute inversion base column-by-column
                 der_temp = 0.
                 height_inv(:,:) = 0
                 iloop3: do i = 1,nx
                        jloop3: do j = 1,ny
                               kloop3: do k = h_cb,nzm
                                         der_temp = (T_diff(i,j,k)/z_diff1(k))
                                         if (der_temp.gt.der_thresh) then
                                             height_inv(i,j) = k
                                             if (der_temp.gt.der_thresh) exit kloop3
                                         else
                                         !for cases where no inv. base was found
                                         height_inv(i,j) = 0
                                         endif
                               enddo kloop3
                               der_temp = 0.
                        enddo jloop3
                 enddo iloop3
                 
                 !------------------------------------------------------------------
                !mcmichael: calculate column-by-column weighted qt and tl for tendency comp.
                !weighted by vertical grid spacing and density

                qt_w(:,:) = 0.
                tl_w(:,:) = 0.
                do i = 1,nx
                        do j = 1,ny
                                !check if inversion base exists
                                if (height_inv(i,j).ne.0) then
                                        rho_tot = SUM(rho(1:height_inv(i,j)))
                                        z_tot = SUM(z_diff1(1:height_inv(i,j)))
                                        do k = 1,height_inv(i,j)
                                                !compute vertical distance and density fraction
                                                !grid correction after z(height_inv(i,j)) is more important for cloud layer avg. 
                                                z_fraction = z_diff1(k) / z_tot
                                                rho_fraction = rho(k) / rho_tot

                                                !calculate grid and density weight
                                                weight = (z_fraction + rho_fraction)*0.5

                                                !accumulate weighted averages
                                                tl_w(i,j) = tl_w(i,j) + t(i,j,k)*weight
                                                qt_w(i,j) = qt_w(i,j) + (qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))*weight
                                        enddo
                                endif
                        enddo
                enddo
         !------------------------------------------------------------------
         !mcmichael: calculate the cloud base and cloud top height
         der_temp = 0.
         iloop1: do i = 1,nx
                jloop1: do j = 1,ny
                       kloop1: do k = 1,nzm
                              der_temp = qcl(i,j,k)
                              if (der_temp.gt.0.) then
                                 height_cb(i,j) = k
                                 if (der_temp.gt.0.) exit kloop1
                              else
                                 !non-cloudy columns set to 1
                                 !to avoid blowups in following section
                                 height_cb(i,j) = 1
                              endif
                       enddo kloop1
                       der_temp = 0.
                enddo jloop1
         enddo iloop1

         der_temp = 0.
         iloop2: do i = 1,nx
                jloop2: do j = 1,ny
                       kloop2: do k = nzm,height_cb(i,j),-1
                              der_temp = qcl(i,j,k)
                              if (der_temp.gt.0.) then
                                 height_ct(i,j) = k
                                 if (der_temp.gt.0.) exit kloop2
                              else
                                 !non-cloudy columns set to 2
                                 !to avoid blowups in following section
                                 height_ct(i,j) = 2
                              endif
                       enddo kloop2
                       der_temp = 0.
                enddo jloop2
         enddo iloop2

         !mcmichael: calculate liquid water path in each column
         !calculate cloud-layer average of q_t and t_l in cloudy regions
         lwp(:,:) = 0.
         tlcl(:,:) = 0.
         qtcl(:,:) = 0.
         do i = 1,nx
                do j = 1,ny
                       !cloud must be at least two grid boxes deep 
                       if (height_ct(i,j).gt.height_cb(i,j)) then
                          rho_tot2 = SUM(rho(height_cb(i,j):height_ct(i,j)))
                          z_tot = SUM(z_diff1(height_cb(i,j):height_ct(i,j)))
                          do k = height_cb(i,j),height_ct(i,j)
                             !compute vertical distance and density fractions
                             z_fraction = z_diff1(k) / z_tot
                             rho_fraction = rho(k) / rho_tot2

                             !calculate grid and density weight
                             weight = (z_fraction + rho_fraction)*0.5

                             lwp(i,j) = lwp(i,j) + rho(k)*qcl(i,j,k)*z_diff1(k)
                             tlcl(i,j) = tlcl(i,j) + t(i,j,k)*weight
                             qtcl(i,j) = qtcl(i,j) + (qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))*weight
                          enddo
                       endif
                enddo
         enddo

         !mcmichael: calculate subcloud layer average of q_t and t_l in cloudy regions
         tlsc(:,:) = 0.
         qtsc(:,:) = 0.
         do i = 1,nx
                do j = 1,ny
                       if (height_cb(i,j).gt.1) then
                          rho_tot3 = SUM(rho(1:height_cb(i,j)))
                          z_tot = SUM(z_diff1(1:height_cb(i,j)))
                          do k = 1,height_cb(i,j)
                             !compute vertical distance and density fractions
                             z_fraction = z_diff1(k) / z_tot   
                             rho_fraction = rho(k) / rho_tot3

                             !calculate grid and density weight
                             weight = (z_fraction + rho_fraction)*0.5

                             tlsc(i,j) = tlsc(i,j) + t(i,j,k)*weight
                             qtsc(i,j) = qtsc(i,j) + (qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))*weight
                          enddo
                       endif
                enddo
         enddo

         !------------------------------------------------------------------

           if(day.lt.shipv2_time0) then
              !if ship has not started injection, don't calculate edge NAc
              NAc_thresh = 0.

           else !ship injection has started
              
                !--------------------------------------------------------------------------
                !after ship injection starts begin calculating NC_env and NAc_thresh 
                
                !grab NC and NAc from microphysics.f90
                Cloudaero(1:nx,1:ny,1:nzm) = Get_cloudaerosol()
                Accumaero(1:nx,1:ny,1:nzm) = Get_accumaerosol()

                !calculate average boundary-layer NAc in each column
                do i = 1,nx
                   do j = 1,ny
                      if (height_inv(i,j).ne.0) then 
                         do k = 1,height_inv(i,j)
                            aero_xy(i,j) = Accumaero(i,j,k)
                            aero_xy_col(i,j) = aero_xy_col(i,j) + aero_xy(i,j)*(z_diff1(k)/z(height_inv(i,j)))
                         enddo
                      endif
                   enddo
                enddo

                !----------------------------------------------------------------------------
                !determine if point is on edge (on edge = 1, in middle = 0)
                !determine x coordinates on current processor to ensure calculations are done
                !on the left and right portions of the domain
                if(dompi) then
                  call task_rank_to_index(rank,it,jt)
                  do i = 1,nx
                     tmpx_pos(i) = dx*(float(i+it))
                  end do
                endif

                edge_factor(:) = 0.
                do i = 1,nx
                   if(tmpx_pos(i).le.((dx*float(nx_gl))*edge_frac).OR.tmpx_pos(i).ge.((dx*float(nx_gl))*(1.0-edge_frac))) then
                     edge_factor(i) = 1.
                   else
                     edge_factor(i) = 0.
                   endif
                enddo

                !point-by-point aerosol averages of edges of domain
                ncenv_numerator = 0.
                ncenv_denominator = 0.
                nacenv_numerator = 0.
                nacenv2_numerator = 0.
                nacenv_denominator = 0.

                do k = 1,nzm
                   rhodz = rho(k)*dz*adz(k)
                   do j = 1,ny
                      do i = 1,nx
                         ncenv_numerator = ncenv_numerator + edge_factor(i)*rhodz*qcl(i,j,k)*cloudaero(i,j,k)
                         ncenv_denominator = ncenv_denominator + edge_factor(i)*rhodz*qcl(i,j,k)
                         if((z(k).lt.z(height_inv(i,j))).and.(height_inv(i,j).ne.0)) then
                            nacenv_numerator = nacenv_numerator + edge_factor(i)*rhodz*accumaero(i,j,k)
                            nacenv2_numerator = nacenv2_numerator + edge_factor(i)*rhodz*accumaero(i,j,k)*accumaero(i,j,k)
                            nacenv_denominator = nacenv_denominator + edge_factor(i)*rhodz
                         end if
                      end do
                   end do
                end do

                !sum across subdomains
                if(dompi) then
                  buffer30(1) = ncenv_numerator
                  buffer30(2) = ncenv_denominator
                  buffer30(3) = nacenv_numerator
                  buffer30(4) = nacenv2_numerator
                  buffer30(5) = nacenv_denominator
                  call task_sum_real(buffer30,buffer31,5) 
                  ncenv_numerator = buffer31(1)
                  ncenv_denominator = buffer31(2)
                  nacenv_numerator = buffer31(3)
                  nacenv2_numerator = buffer31(4)
                  nacenv_denominator = buffer31(5)
                end if

                NC_env = ncenv_numerator / (1.e-16 + ncenv_denominator)
                NAC_mean2 = nacenv_numerator / (1.e-16 + nacenv_denominator)
                NAC2_mean = nacenv2_numerator / (1.e-16 + nacenv_denominator)
                NAC_std2 = sqrt(NAC2_mean - NAC_mean2**2)
               
                NAc_thresh = NAC_mean2 + 3.0*NAC_std2

                if(masterproc) then
                   print*, 'NC in environment (#/mg) =', NC_env*1.e-6
                   print*, 'NAc in environment (#/mg) =', NAC_mean2*1.e-6
                   print*, 'NAc standard deviation (#/mg) =', NAC_std2*1.e-6
                endif            

           endif
         endif
         !-----------------------------------------------------------------------------------------------------------------

         condavg_mask(:,:,:,:) = 0.

         do k = 1,nzm
            if(LES) then
               coef=0.
            else
               coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
            endif
            do j = 1,ny
               do i = 1,nx
                  
                  if((icondavg_cld.gt.0).and.(qcl(i,j,k)+qci(i,j,k).gt.coef)) then
                     condavg_mask(i,j,k,icondavg_cld) = 1. ! cloud
                  end if

                  if(icondavg_cor.gt.0) then
                     ! updraft (w>1) core (tv'>0) statistics
                     ! in LES, buoyant cloudy statistics
                     condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef
                     condition = tvirt(i,j,k).gt.tvz(k) 
                     if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).gt.2.
                     if(LES) condition=condition_cl.and.condition
                     if(condition) condavg_mask(i,j,k,icondavg_cor) = 1. ! core
                  end if

                  if(icondavg_cordn.gt.0) then
                     ! downdraft (w<-1) core (tv'>0) statistics
                     ! in LES, buoyant, saturated or rainy statistics
                     condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef &
                                .or. qpl(i,j,k)+qpi(i,j,k).gt.1.e-4 
                     condition = tvirt(i,j,k).lt.tvz(k) 
                     if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).lt.-2.
                     if(LES) condition=condition_cl.and.condition
                     if(condition) condavg_mask(i,j,k,icondavg_cordn) = 1. ! downdraft core
                  end if


                  condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef
                  if((icondavg_satup.gt.0).AND. &
                       (condition_cl.AND.w(i,j,k)+w(i,j,k+1).ge.0.)) then
                     condavg_mask(i,j,k,icondavg_satup) = 1. ! saturated updraft
                  end if

                  if((icondavg_satdn.gt.0).AND. &
                       (condition_cl.AND.w(i,j,k)+w(i,j,k+1).lt.0.)) then
                     condavg_mask(i,j,k,icondavg_satdn) = 1. ! saturated downdraft
                  end if

                  if((icondavg_env.gt.0).AND.(.NOT.condition_cl)) then
                     condavg_mask(i,j,k,icondavg_env) = 1. ! cloud-free environment
                  end if
                  
		  !mcmichael: only loop through i,j indices for ship track stats
                  if(doShipTrackConditionals) then
                     if(k.eq.1) then
		  
                  	!conditional averages for ship track 
                  	!computed for mixed-layer model profiles
                  	!make entire column 1 if cloud is present
			qclz(:) = qcl(i,j,:)
			qciz(:) = qci(i,j,:)
                  	tmpqcl = MAXVAL(qclz(:))
                  	tmpqci = MAXVAL(qciz(:))
			!if any liquid is present = .true.
                  	condition_cloud = (tmpqcl + tmpqci).gt.0
                  
                  	!Calculate ship track condition                    
                  	if(day.lt.shipv2_time0) then
                     		condition_track = .false.
                  	else
                     		condition_track = aero_xy_col(i,j).gt.NAc_thresh.and.NAc_thresh.ne.0.
                  	endif

                  	!Print out aerosol threshold to log file
                  	if(masterproc.and.(i+j+k.eq.3)) then 
                    		print*, 'Aerosol threshold for ship track (#/mg) =', NAc_thresh*1.e-6
                  	endif
                  
                        !added check for inversion being present 
                  	if((icondavg_sh_cloud.gt.0).AND. &
                       		(condition_cloud.AND.condition_track) &
                                .AND.(height_inv(i,j).ne.0)) then
                     		condavg_mask(i,j,:,icondavg_sh_cloud) = 1. !ship and cloud
                  	end if
                 
                  	if((icondavg_sh_clear.gt.0).AND.(.NOT.condition_cloud).AND. &
                       		(condition_track) &
                                .AND.(height_inv(i,j).ne.0)) then
                     		condavg_mask(i,j,:,icondavg_sh_clear) = 1. !ship and clear
                  	end if

                  	if((icondavg_no_sh_cloud.gt.0).AND. &
                       		(condition_cloud).AND.(.NOT.condition_track) &
                                .AND.(height_inv(i,j).ne.0)) then
                     		condavg_mask(i,j,:,icondavg_no_sh_cloud) = 1. !no ship and cloud
                  	end if                  

                  	if((icondavg_no_sh_clear.gt.0).AND. &
                       		(.NOT.condition_cloud).AND.(.NOT.condition_track) &
                                .AND.(height_inv(i,j).ne.0)) then
                     		condavg_mask(i,j,:,icondavg_no_sh_clear) = 1. !no ship, no cloud
                  	end if
			
		    end if
                  end if

               end do
            end do
         end do

         do ncond = 1,ncondavg
            cld(:) = 0.
            wcl(:) = 0.
            ucl(:) = 0.
            vcl(:) = 0.
            wacl(:) = 0.
            tcl(:) = 0.
            tacl(:) = 0.
            tvcl(:)= 0.
            tvcla(:)= 0.
            qcll(:) = 0.
            qccl(:)= 0.
            qicl(:)= 0.
            qpcl(:)= 0.
            tvwcl(:)= 0.
            twcl(:)= 0.
            qwcl(:)= 0.
            qcwcl(:)= 0.
            qiwcl(:)= 0.
            dse(:)=0.
            mse(:)=0.
            sse(:)=0.
 
            !mcmichael: mass budget/MLM budget terms
            zb(:) = 0.
            zbct(:) = 0.
            wlo(:) = 0.
            wla(:) = 0.
            wloct(:) = 0.
            wlact(:) = 0.
            swu0(:) = 0.
            swd0(:) = 0.
            lwu0(:) = 0.
            lwd0(:) = 0.
            swut(:) = 0.
            swdt(:) = 0.
            lwut(:) = 0.
            lwdt(:) = 0.
            rdiv(:) = 0.
            ldiv(:) = 0.
            sdiv(:) = 0.
            swutct(:) = 0.
            swdtct(:) = 0.
            lwutct(:) = 0.
            lwdtct(:) = 0.
            rdivct(:) = 0.
            ldivct(:) = 0.
            sdivct(:) = 0.
            sh_f(:) = 0.
            lh_f(:) = 0.
            pr_s(:) = 0.
            pr_z(:) = 0.
            pr_zct(:) = 0.
            pr_zb(:) = 0.
            tl_avg(:) = 0.
            tlcl_avg(:) = 0.
            qtcl_avg(:) = 0.
            tlsc_avg(:) = 0.
            qtsc_avg(:) = 0.
            tcb(:) = 0.
            qt_avg(:) = 0.
            cb_h_avg(:) = 0.
            cl_d_avg(:) = 0.
            lwp_avg(:) = 0.

            !mcmichael: advection-related variables
            u_adv_t(:) = 0.
            u_adv_q(:) = 0.
            v_adv_t(:) = 0.
            v_adv_q(:) = 0.
            w_adv_t(:) = 0.
            w_adv_q(:) = 0.

            !mcmichael: radiation streams sampled conditionally
            sw_u(:) = 0.
            sw_d(:) = 0.
            lw_u(:) = 0.
            lw_d(:) = 0.

            !mcmichael: variances for ship track stats
            ts2(:) = 0.
            qs2(:) = 0.
            ws2(:) = 0.
            us2(:) = 0.
            vs2(:) = 0.

            !mcmichael: vertical turbulent fluxes for ship tracks
            wtl(:) = 0.
            wqt(:) = 0.
            wtv(:) = 0.

            !mcmichael: precipitation flux for ship tracks
            pfs(:) = 0.

            !mcmichael: relaxation timescale variables
            disstke(:) = 0.
            difftke(:) = 0.
            restke(:) = 0.

            !bloss: conditional u,v anomalies
            ucla(:) = 0. 
            vcla(:) = 0.

            !bloss: pressure gradients
            dpdxcl(:) = 0.
            dpdycl(:) = 0.
            dpdzcl(:) = 0.

            !bloss: add momentum fluxes
            uwsbcl(:) = 0. 
            vwsbcl(:) = 0.
            uwlecl(:) = 0.
            vwlecl(:) = 0.

            !bloss: frozen moist static energy
            fmse(:) = 0.
            fmsecla(:) = 0.

            !bloss: mass flux and mass-flux weighted stats in conditional category
            rhowcl(:) = 0.

            rhowmsecl(:) = 0.
            rhowtlcl(:) = 0. 
            rhowqtcl(:) = 0. 
            rhowtvcl(:) = 0. 

            rhowmsecla(:) = 0.
            rhowqtcla(:) = 0. 
            rhowtlcla(:) = 0. 
            rhowtvcla(:) = 0. 

            rhouwcl(:) = 0.
            rhovwcl(:) = 0.
            rhowwcl(:) = 0.

            do k=1,nzm
               if(LES) then
                  coef=0.
               else
                  coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
               endif
               kb = max(1,k-1)
               kc = min(nzm,k+1)
               do j=1,ny
                  jb = YES3D*(j-1) + (1-YES3D)
                  do i=1,nx

                     if(condavg_mask(i,j,k,ncond).gt.0) then

                        ! gather conditional statistics
                        cld(k)=cld(k) + 1
                        tmp(1)=0.5*(w(i,j,k+1)+w(i,j,k))
                        wcl(k) = wcl(k) + tmp(1)
                        ucl(k) = ucl(k) + u(i,j,k) + ug !bloss: include ground speed
                        vcl(k) = vcl(k) + v(i,j,k) + vg                        

                        if(doShipTrackConditionals.and.ncond.ge.7) then !mcmichael 
                                
                           if(k.eq.1) then
                              !cloud top (ct) statistics are only valid for
                              !cloudy conditionals ST1, ST3

                              if (height_inv(i,j).ne.0) then
                                 zb(:) = zb(:) + z(height_inv(i,j)) !inversion base height
                                 wlo(:) = wlo(:) + w(i,j,height_inv(i,j)) !regional w
                                 wla(:) = wla(:) + wsub(height_inv(i,j)) !large-scale w
                                 swut(:) = swut(:) + swUp3D(i,j,height_inv(i,j))
                                 swdt(:) = swdt(:) + swDown3D(i,j,height_inv(i,j))
                                 swu0(:) = swu0(:) + swUp3D(i,j,1)
                                 swd0(:) = swd0(:) + swDown3D(i,j,1)
                                 lwu0(:) = lwu0(:) + lwUp3D(i,j,1)
                                 lwd0(:) = lwd0(:) + lwDown3D(i,j,1)
                                 lwut(:) = lwut(:) + lwUp3D(i,j,height_inv(i,j))
                                 lwdt(:) = lwdt(:) + lwDown3D(i,j,height_inv(i,j))
                                 ldiv(:) = (lwdt(:) - lwut(:)) - (lwd0(:) - lwu0(:))
                                 sdiv(:) = (swdt(:) - swut(:)) - (swd0(:) - swu0(:))
                                 rdiv(:) = (swdt(:) - swut(:) + lwdt(:) - lwut(:)) &
                                        -(swd0(:) - swu0(:) + lwd0(:) - lwu0(:))
                                 pr_z(:) = pr_z(:) + 0.5*(w(i,j,height_inv(i,j)+1) &
                                        +w(i,j,height_inv(i,j)))*(qpl(i,j,height_inv(i,j)) &
                                        -qrz(height_inv(i,j))) !inv. top prec. flux
                              endif

                              zbct(:) = zbct(:) + z(height_ct(i,j)) !cloud-top height                  
                              wloct(:) = wloct(:) + w(i,j,height_ct(i,j)) !regional w
                              wlact(:) = wlact(:) + wsub(height_ct(i,j)) !large-scale w
                              swutct(:) = swutct(:) + swUp3D(i,j,height_ct(i,j))
                              swdtct(:) = swdtct(:) + swDown3D(i,j,height_ct(i,j))
                              lwutct(:) = lwutct(:) + lwUp3D(i,j,height_ct(i,j))
                              lwdtct(:) = lwdtct(:) + lwDown3D(i,j,height_ct(i,j))
                              rdivct(:) = (swdtct(:) - swutct(:) + lwdtct(:) - lwutct(:)) &
                                          -(swd0(:) - swu0(:) + lwd0(:) - lwu0(:))
                              ldivct(:) = (lwdtct(:) - lwutct(:)) - (lwd0(:) - lwu0(:))
                              sdivct(:) = (swdtct(:) - swutct(:)) - (swd0(:) - swu0(:))
                              sh_f(:) = sh_f(:) + fluxbt(i,j) !SHF at sfc
                              lh_f(:) = lh_f(:) + fluxbq(i,j) !LHF at sfc
                              pr_s(:) = pr_s(:) + 0.5*(w(i,j,k+1)+w(i,j,k)) &
                                        *(qpl(i,j,1)-qrz(1)) !sfc prec. flux
                              pr_zct(:) = pr_zct(:) + 0.5*(w(i,j,height_ct(i,j)+1) &
                                          +w(i,j,height_ct(i,j)))*(qpl(i,j,height_ct(i,j)) &
                                          -qrz(height_ct(i,j))) !inv. top prec. flux
                              pr_zb(:) = pr_zb(:) + 0.5*(w(i,j,height_cb(i,j)+1) &
                                          +w(i,j,height_cb(i,j)))*(qpl(i,j,height_cb(i,j)) &
                                          -qrz(height_cb(i,j))) !cloud base prec. flux
                              tl_avg(:) = tl_avg(:) + tl_w(i,j) !weighted tl average
                              tcb(:) = tcb(:) + tabs(i,j,height_cb(i,j)) !cb temp
                              qt_avg(:) = qt_avg(:) + qt_w(i,j) !weighted qt average
                              cb_h_avg(:) = cb_h_avg(:) + z(height_cb(i,j)) 
                              cl_d_avg(:) = cl_d_avg(:) + z(height_ct(i,j)) - z(height_cb(i,j))
                              lwp_avg(:) = lwp_avg(:) + lwp(i,j) 
                              tlcl_avg(:) = tlcl_avg(:) + tlcl(i,j)
                              tlsc_avg(:) = tlsc_avg(:) + tlsc(i,j)
                              qtcl_avg(:) = qtcl_avg(:) + qtcl(i,j)
                              qtsc_avg(:) = qtsc_avg(:) + qtsc(i,j)
                           endif     

                           sw_u(k) = sw_u(k) + swUp3D(i,j,k) !upwelling sw
                           sw_d(k) = sw_d(k) + swDown3D(i,j,k) !downwelling sw
                           lw_u(k) = lw_u(k) + lwUp3D(i,j,k) !upwelling lw
                           lw_d(k) = lw_d(k) + lwDown3D(i,j,k) !downwelling lw
                        
                           ts2(k) = ts2(k) + (t(i,j,k)-t0(k))**2 !tl variance
                           qs2(k) = qs2(k) + (qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k))**2 !qt variance
                           ws2(k) = ws2(k) + 0.5*(w(i,j,k+1)**2+w(i,j,k)**2) !w variance 
                           us2(k) = us2(k) + (u(i,j,k) - u0(k))**2 !u variance
                           vs2(k) = vs2(k) + (v(i,j,k) - v0(k))**2 !v variance

                           pfs(k) = pfs(k)+tmp(1)*(qpl(i,j,k)-qrz(k)) !liquid precipitation flux
     
                           wtl(k) = wtl(k) + 0.5*w(i,j,k)*(t(i,j,kb)-t0(kb)+t(i,j,k)-t0(k)) !TL flux   
                           wqt(k) = wqt(k) + 0.5*w(i,j,k)*(qv(i,j,kb)+qcl(i,j,kb)+qci(i,j,kb)-q0(kb)+ &
                                    qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k)) !QT flux
                           wtv(k) = wtv(k) + 0.5*w(i,j,k)*(tvirt(i,j,kb)-tvz(kb)+tvirt(i,j,k)-tvz(k)) !TV flux

                           !compute full advection profile
                           u_adv_t(k) = u_adv_t(k) + ( u(i,j,k)*(0.5*(t(i,j,k) + t(i-1,j,k))) &
                                                     - u(i+1,j,k)*(0.5*(t(i+1,j,k) + t(i,j,k))))/dx
                           v_adv_t(k) = v_adv_t(k) + ( v(i,j,k)*(0.5*(t(i,j,k) + t(i,j-1,k))) &
                                                     - v(i,j+1,k)*(0.5*(t(i,j+1,k) + t(i,j,k))))/dy
                           u_adv_q(k) = u_adv_q(k) + ( u(i,j,k)*(0.5*(tot_w(i,j,k) + tot_w(i-1,j,k))) &
                                                     - u(i+1,j,k)*(0.5*(tot_w(i+1,j,k) + tot_w(i,j,k))))/dx
                           v_adv_q(k) = v_adv_q(k) + ( v(i,j,k)*(0.5*(tot_w(i,j,k) + tot_w(i,j-1,k))) &
                                                     - v(i,j+1,k)*(0.5*(tot_w(i,j+1,k) + tot_w(i,j,k))))/dy
                           w_adv_t(k) = w_adv_t(k) + ( rhow(k)*w(i,j,k)*(0.5*(t(i,j,k) + t(i,j,kb))) &
                                                     - rhow(kc)*w(i,j,kc)*(0.5*(t(i,j,kc) + t(i,j,k))))/(rho(k)*z_diff1(k))
                           w_adv_q(k) = w_adv_q(k) + ( rhow(k)*w(i,j,k)*(0.5*(tot_w(i,j,k) + tot_w(i,j,kb))) &
                                                     - rhow(kc)*w(i,j,kc)*(0.5*(tot_w(i,j,kc) + tot_w(i,j,k))))/(rho(k)*z_diff1(k))
						   

                           !relaxation timescale variables
                           disstke(k) = disstke(k) + tkediss(i,j,k)
                           difftke(k) = difftke(k) + tkediff(i,j,k)
                           restke(k) = restke(k) + 0.5*(((u(i,j,k) - u0(k))**2) + ((v(i,j,k) - v0(k))**2) &
                                       + 0.5*(w(i,j,k+1)**2+w(i,j,k)**2))

                        endif
 
                        ucla(k) = ucla(k) + u(i,j,k) - u0(k) !bloss: u,v anomalies
                        vcla(k) = vcla(k) + v(i,j,k) - v0(k)
                        qcc=qcl(i,j,k)
                        qii=qci(i,j,k)
                        dse(k)=dse(k)+tabs(i,j,k)+gamaz(k)	
                        mse(k)=mse(k)+tabs(i,j,k)+gamaz(k)+fac_cond*qv(i,j,k)	
                        tcl(k) = tcl(k) + t(i,j,k)
                        qcll(k) = qcll(k) + (qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))
                        qccl(k) = qccl(k) + qcc
                        qicl(k) = qicl(k) + qii
                        qpcl(k) = qpcl(k) + qpl(i,j,k) + qpi(i,j,k)
                        tvcl(k) = tvcl(k) + tvirt(i,j,k)	 
                        tvcla(k) = tvcla(k) + tvirt(i,j,k) - tvz(k)	 
                        tacl(k) = tacl(k) + tabs(i,j,k)
                        twcl(k) = twcl(k) + t(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
                        qwcl(k) = qwcl(k) + (qv(i,j,k)+qcl(i,j,k)+qci(i,j,k))*0.5*(w(i,j,k+1)+w(i,j,k))
                        tvwcl(k) = tvwcl(k)+tvirt(i,j,k)*0.5*(w(i,j,k+1)+w(i,j,k))
                        qcwcl(k) = qcwcl(k) + qcc*0.5*(w(i,j,k+1)+w(i,j,k))
                        qiwcl(k) = qiwcl(k) + qii*0.5*(w(i,j,k+1)+w(i,j,k))

                        !bloss: frozen MSE and anomaly
                        fmse(k)=fmse(k)+t(i,j,k)+fac_cond*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)+qpl(i,j,k)+qpi(i,j,k))
                        fmsecla(k)=fmsecla(k)+t(i,j,k)-t0(k) &
                             +fac_cond*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)+qpl(i,j,k)+qpi(i,j,k)-q0(k)-qp0(k))

                        !bloss: pressure gradient forces
                        dpdxcl(k) = dpdxcl(k) - (p(i,j,k)-p(i-1,j,k))/(dx*rho(k))
                        dpdycl(k) = dpdycl(k) - (p(i,j,k)-p(i,jb,k))/(dy*rho(k))
                        dpdzcl(k) = dpdzcl(k) &
                             - 0.5*(p(i,j,kc)/rho(kc)-p(i,j,k)/rho(k))/(dz*adzw(kc)) &
                             - 0.5*(p(i,j,k)/rho(k)-p(i,j,kb)/rho(kb))/(dz*adzw(k))

                        !bloss: add momentum fluxes
                        if(k.eq.1) then
                           ! surface momentum flux anomaly
                           uwsubgrid = fluxbu(i,j) !surface momentum flux/drag
                           vwsubgrid = fluxbv(i,j) !surface momentum flux/drag

                           uwresolved = 0. ! no resolved momentum flux at surface
                           vwresolved = 0. ! no resolved momentum flux at surface

                        else              
                           ! momentum flux anomaly above surface

                           ! subgrid
               !            uwsubgrid = -(0.25*grdf_z(k-1)/dz) &
               !                 *(tk(i,j,k-1)+tk(i-1,j,k-1)+tk(i,j,k)+tk(i-1,j,k)) &
               !                 *( (u(i,j,k)-u(i,j,k-1))/adzw(k) &
               !                   + (w(i,j,k)-w(i-1,j,k))*dz/dx)
               !            vwsubgrid = -(0.25*grdf_z(k-1)/dz) &
               !                 *(tk(i,j,k-1)+tk(i,jb,k-1)+tk(i,j,k)+tk(i,jb,k)) &
               !                 *( (v(i,j,k)-v(i,j,k-1))/adzw(k) &
               !                   + (w(i,j,k)-w(i,jb,k))*dz/dy)

                            uwsubgrid = 0.
                            vwsubgrid = 0.
                           ! resolved
                           uwresolved = 0.25*(w(i,j,k)+w(i-1,j,k)) &
                                *(u(i,j,k)+u(i,j,k-1)-u0(k)-u0(k-1))
                           vwresolved = 0.25*(w(i,j,k)+w(i,jb,k)) &
                                *(v(i,j,k)+v(i,j,k-1)-v0(k)-v0(k-1))
                        end if

                        uwsbcl(k) = uwsbcl(k) + uwsubgrid
                        vwsbcl(k) = vwsbcl(k) + vwsubgrid
                        uwlecl(k) = uwlecl(k) + uwresolved + uwsubgrid
                        vwlecl(k) = vwlecl(k) + vwresolved + vwsubgrid

                        !bloss: add mass flux and mass flux weighted stats
                        tmprhow = 0.5*rho(k)*(w(i,j,k+1) + w(i,j,k))
                        tmpmse = t(i,j,k) &
                             + fac_cond*(qv(i,j,k) + qcl(i,j,k) + qci(i,j,k) &
                             + qpl(i,j,k) + qpi(i,j,k))
                        tmpqt = qv(i,j,k) + qcl(i,j,k) + qci(i,j,k)

                        rhowcl(k) = rhowcl(k) + tmprhow
                        rhowmsecl(k) = rhowmsecl(k) + tmprhow*tmpmse
                        rhowmsecla(k) = rhowmsecla(k) &
                             + tmprhow*(tmpmse - t0(k) - fac_cond*(q0(k) + qp0(k)))
                        rhowqtcl(k) = rhowqtcl(k) + tmprhow*tmpqt
                        rhowqtcla(k) = rhowqtcla(k) + tmprhow*(tmpqt-q0(k))
                        rhowtlcl(k) = rhowtlcl(k) + tmprhow*(t(i,j,k))
                        rhowtlcla(k) = rhowtlcla(k) + tmprhow*(t(i,j,k)-t0(k))
                        rhowtvcl(k) = rhowtvcl(k) + tmprhow*tvirt(i,j,k)
                        rhowtvcla(k) = rhowtvcla(k) + tmprhow*(tvirt(i,j,k)-tvz(k))

                        rhouwcl(k) = rhouwcl(k) + tmprhow*(u(i,j,k) - u0(k))
                        rhovwcl(k) = rhovwcl(k) + tmprhow*(v(i,j,k) - v0(k))
                        rhowwcl(k) = rhowwcl(k) + rho(k)*0.5*(w(i,j,k)**2 + w(i,j,k+1)**2)

                     endif
                  end do                  
               end do
               condavg_factor(k,ncond) = condavg_factor(k,ncond)+cld(k)
               wacl(k) = wcl(k)

               if(doShipTrackConditionals.and.ncond.ge.7) then
                  wtl(k) = wtl(k)*rho(k)*cp
                  wqt(k) = wqt(k)*rho(k)*lcond
                  wtv(k) = wtv(k)*rho(k)*cp
                  sh_f(k) = sh_f(k)*rhow(k)*cp
                  lh_f(k) = lh_f(k)*rhow(k)*lcond
               endif

            end do

            call hbuf_put(TRIM(condavgname(ncond)),cld,factor_xy)
            call hbuf_put('W'//TRIM(condavgname(ncond)),wcl,1.)
            call hbuf_put('U'//TRIM(condavgname(ncond)),ucl,1.)
            call hbuf_put('V'//TRIM(condavgname(ncond)),vcl,1.)

            if(doShipTrackConditionals.and.ncond.ge.7) then

                !ship track conditionals for mixed-layer model
                call hbuf_put('ZB'//TRIM(condavgname(ncond)),zb,1.)
                call hbuf_put('ZBCT'//TRIM(condavgname(ncond)),zbct,1.) 
                call hbuf_put('WLO'//TRIM(condavgname(ncond)),wlo,1.)
                call hbuf_put('WLA'//TRIM(condavgname(ncond)),wla,1.)
                call hbuf_put('WLOT'//TRIM(condavgname(ncond)),wloct,1.)
                call hbuf_put('WLAT'//TRIM(condavgname(ncond)),wlact,1.)
                call hbuf_put('RDIV'//TRIM(condavgname(ncond)),rdiv,1.)
                call hbuf_put('LDIV'//TRIM(condavgname(ncond)),ldiv,1.)
                call hbuf_put('SDIV'//TRIM(condavgname(ncond)),sdiv,1.)
                call hbuf_put('RDVT'//TRIM(condavgname(ncond)),rdivct,1.)
                call hbuf_put('LDVT'//TRIM(condavgname(ncond)),ldivct,1.)
                call hbuf_put('SDVT'//TRIM(condavgname(ncond)),sdivct,1.)
                call hbuf_put('SH_F'//TRIM(condavgname(ncond)),sh_f,1.)
                call hbuf_put('LH_F'//TRIM(condavgname(ncond)),lh_f,1.)
                call hbuf_put('PR_S'//TRIM(condavgname(ncond)),pr_s,1.)
                call hbuf_put('PR_Z'//TRIM(condavgname(ncond)),pr_z,1.)
                call hbuf_put('PRZT'//TRIM(condavgname(ncond)),pr_zct,1.)
                call hbuf_put('PRZB'//TRIM(condavgname(ncond)),pr_zb,1.)
                call hbuf_put('QT_A'//TRIM(condavgname(ncond)),qt_avg,1.)
                call hbuf_put('TL_A'//TRIM(condavgname(ncond)),tl_avg,1.)
                call hbuf_put('LWP'//TRIM(condavgname(ncond)),lwp_avg,1.)
                call hbuf_put('TCB'//TRIM(condavgname(ncond)),tcb,1.)
                call hbuf_put('CBH'//TRIM(condavgname(ncond)),cb_h_avg,1.)
                call hbuf_put('TLCL'//TRIM(condavgname(ncond)),tlcl_avg,1.)
                call hbuf_put('TLSC'//TRIM(condavgname(ncond)),tlsc_avg,1.)
                call hbuf_put('QTCL'//TRIM(condavgname(ncond)),qtcl_avg,1.)
                call hbuf_put('QTSC'//TRIM(condavgname(ncond)),qtsc_avg,1.)
                call hbuf_put('SW_U'//TRIM(condavgname(ncond)),sw_u,1.)
                call hbuf_put('SW_D'//TRIM(condavgname(ncond)),sw_d,1.)
                call hbuf_put('LW_U'//TRIM(condavgname(ncond)),lw_u,1.)
                call hbuf_put('LW_D'//TRIM(condavgname(ncond)),lw_d,1.)
                call hbuf_put('TS2'//TRIM(condavgname(ncond)),ts2,1.)
                call hbuf_put('QS2'//TRIM(condavgname(ncond)),qs2,1.e6)
                call hbuf_put('WS2'//TRIM(condavgname(ncond)),ws2,1.)
                call hbuf_put('US2'//TRIM(condavgname(ncond)),us2,1.)
                call hbuf_put('VS2'//TRIM(condavgname(ncond)),vs2,1.)
                call hbuf_put('WTL'//TRIM(condavgname(ncond)),wtl,1.)
                call hbuf_put('WQT'//TRIM(condavgname(ncond)),wqt,1.)
                call hbuf_put('WTV'//TRIM(condavgname(ncond)),wtv,1.) 
                call hbuf_put('PFS'//TRIM(condavgname(ncond)),pfs,1.)
                call hbuf_put('UADT'//TRIM(condavgname(ncond)),u_adv_t,1.)
                call hbuf_put('UADQ'//TRIM(condavgname(ncond)),u_adv_q,1.)
                call hbuf_put('VADT'//TRIM(condavgname(ncond)),v_adv_t,1.)
                call hbuf_put('VADQ'//TRIM(condavgname(ncond)),v_adv_q,1.)
                call hbuf_put('WADT'//TRIM(condavgname(ncond)),w_adv_t,1.)
                call hbuf_put('WADQ'//TRIM(condavgname(ncond)),w_adv_q,1.)
                call hbuf_put('DIS'//TRIM(condavgname(ncond)),disstke,1.)
                call hbuf_put('DIF'//TRIM(condavgname(ncond)),difftke,1.)
                call hbuf_put('TKER'//TRIM(condavgname(ncond)),restke,1.)
                !---------------------------------------------
                 
            endif

            call hbuf_put('DSE'//TRIM(condavgname(ncond)),dse,1.)
            call hbuf_put('MSE'//TRIM(condavgname(ncond)),mse,1.)
            call hbuf_put('TL'//TRIM(condavgname(ncond)),tcl,1.)
            call hbuf_put('TV'//TRIM(condavgname(ncond)),tvcl,1.)
            call hbuf_put('TV'//TRIM(condavgname(ncond))//'A',tvcla,1.)
            call hbuf_put('TA'//TRIM(condavgname(ncond)),tacl,1.)
            call hbuf_put('QT'//TRIM(condavgname(ncond)),qcll,1.e3)
            !bloss               call hbuf_put('QC'//TRIM(condavgname(ncond)),qccl,1.e3)
            !bloss               call hbuf_put('QI'//TRIM(condavgname(ncond)),qicl,1.e3)
            call hbuf_put('QN'//TRIM(condavgname(ncond)),qccl+qicl,1.e3)
            call hbuf_put('QP'//TRIM(condavgname(ncond)),qpcl,1.e3)
            call hbuf_put('W'//TRIM(condavgname(ncond))//'A',wacl,factor_xy)
            call hbuf_put('TLW'//TRIM(condavgname(ncond)),twcl,factor_xy)
            call hbuf_put('TVW'//TRIM(condavgname(ncond)),tvwcl,factor_xy)
            call hbuf_put('QTW'//TRIM(condavgname(ncond)),qwcl,factor_xy*1.e3)
            call hbuf_put('QCW'//TRIM(condavgname(ncond)),qcwcl,factor_xy*1.e3)
            call hbuf_put('QIW'//TRIM(condavgname(ncond)),qiwcl,factor_xy*1.e3)

            !bloss: add mass flux and mass-flux weighted MSE/QT/TV and anomalies
            call hbuf_put('MF'//TRIM(condavgname(ncond)),rhowcl,factor_xy)

            call hbuf_put('MFH'//TRIM(condavgname(ncond)),rhowmsecl,factor_xy)
            call hbuf_put('MFTL'//TRIM(condavgname(ncond)),rhowtlcl,factor_xy)
            call hbuf_put('MFQT'//TRIM(condavgname(ncond)),rhowqtcl,factor_xy*1.e3)
            call hbuf_put('MFTV'//TRIM(condavgname(ncond)),rhowtvcl,factor_xy)

            call hbuf_put('RUW'//TRIM(condavgname(ncond)),rhouwcl,factor_xy)
            call hbuf_put('RVW'//TRIM(condavgname(ncond)),rhovwcl,factor_xy)
            call hbuf_put('RWW'//TRIM(condavgname(ncond)),rhowwcl,factor_xy)

            call hbuf_put('MFH'//TRIM(condavgname(ncond))//'A',rhowmsecla,factor_xy)
            call hbuf_put('MFTL'//TRIM(condavgname(ncond))//'A',rhowtlcla,factor_xy)
            call hbuf_put('MFQT'//TRIM(condavgname(ncond))//'A',rhowqtcla,factor_xy*1.e3)
            call hbuf_put('MFTV'//TRIM(condavgname(ncond))//'A',rhowtvcla,factor_xy)

            !bloss: add momentum fluxes and horizontal velocity anomalies
            call hbuf_put('UW'//TRIM(condavgname(ncond)),uwlecl,1.)
            call hbuf_put('VW'//TRIM(condavgname(ncond)),vwlecl,1.)
            call hbuf_put('UWSB'//TRIM(condavgname(ncond)),uwsbcl,1.)
            call hbuf_put('VWSB'//TRIM(condavgname(ncond)),vwsbcl,1.)
            call hbuf_put('U'//TRIM(condavgname(ncond))//'A',ucla,1.)
            call hbuf_put('V'//TRIM(condavgname(ncond))//'A',vcla,1.)

            !bloss: frozen moist static energy
            call hbuf_put('HF'//TRIM(condavgname(ncond)),fmse,1.)
            call hbuf_put('HF'//TRIM(condavgname(ncond))//'A',fmsecla,1.)

            !bloss: pressure gradient forces
            call hbuf_put('UPGF'//TRIM(condavgname(ncond)),dpdxcl,1.)
            call hbuf_put('VPGF'//TRIM(condavgname(ncond)),dpdycl,1.)
            call hbuf_put('WPGF'//TRIM(condavgname(ncond)),dpdzcl,1.)

         end do ! ncond = 1,ncondstats

!-------------------------------------------------------------
!	Mass flux, hydrometeor fraction statistics
!-------------------------------------------------------------

	do k=1,nzm

	 hydro(k) = 0.
	 prof1(k)=0.
	 prof2(k)=0.
	 prof3(k)=0.
	 prof4(k)=0.
	 if(LES) then
	  coef=0.
	 else
	  coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
	 endif
	 do j=1,ny
	  do i=1,nx
	    if(qcl(i,j,k)+qci(i,j,k).gt.coef) then
	      hydro(k) = hydro(k) + 1
	      tmp(1)=0.5*(w(i,j,k+1)+w(i,j,k))
	      if(tmp(1).gt.0.) then
		prof1(k)=prof1(k)+rho(k)*tmp(1)
	      else
	        prof2(k)=prof2(k)+rho(k)*tmp(1)
	      endif	
	    elseif(qpl(i,j,k)+qpi(i,j,k).gt.1.e-4) then
	      hydro(k) = hydro(k) + 1
	      if(w(i,j,k)+w(i,j,k+1).lt.0.) &
 	         prof3(k)=prof3(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	    endif
	  end do
	 end do
	 prof4(k)=prof1(k)+prof2(k)+prof3(k)

	end do
		
	call hbuf_put('HYDRO',hydro,factor_xy)
	call hbuf_put('MCUP',prof1,factor_xy)
	call hbuf_put('MCDNS',prof2,factor_xy)
	call hbuf_put('MCDNU',prof3,factor_xy)
	call hbuf_put('MC',prof4,factor_xy)

!-------------------------------------------------------------
!	Updraft Core statistics:
!-------------------------------------------------------------


	do k=1,nzm
	 cldd(k) = 0.
	 prof1(k)=0.
         if(LES) then
          coef=0.
         else
          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
         endif
	 do j=1,ny
	  do i=1,nx
	    condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef
     	    condition = tvirt(i,j,k).gt.tvz(k) 
	    if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).gt.2.
            if(LES) condition=condition_cl.and.condition
	    if(condition) then
	      prof1(k)=prof1(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	      if(condition_cl) then
	        cldd(k)=cldd(k)+1
	      end if
	    endif
	  end do
	 end do
	end do
		
	call hbuf_put('CORECL',cldd,factor_xy)
!-------------------------------------------------------------
!	Cloud Downdraft Core statistics:
!-------------------------------------------------------------


	do k=1,nzm
	 cldd(k) = 0.
	 prof2(k)=0.
	 prof3(k)=0.
         if(LES) then
          coef=0.
         else
          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
         endif
	 do j=1,ny
	  do i=1,nx
	    condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef .or. qpl(i,j,k)+qpi(i,j,k).gt.1.e-4 
     	    condition = tvirt(i,j,k).lt.tvz(k) 
	    if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).lt.-2.
	    if(LES) condition=condition_cl.and.condition
	    if(condition) then
	      if(condition_cl) then
	        prof2(k)=prof2(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	        cldd(k)=cldd(k) + 1
	      else
	        prof3(k)=prof3(k)+rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))      
	      end if
	    endif
	  end do
	 end do
	 prof4(k)=prof1(k)+prof2(k)+prof3(k)

	end do
		
	call hbuf_put('COREDNCL',cldd,factor_xy)
	call hbuf_put('MCRUP',prof1,factor_xy)
	call hbuf_put('MCRDNS',prof2,factor_xy)
	call hbuf_put('MCRDNU',prof3,factor_xy)
	call hbuf_put('MCR',prof4,factor_xy)

!---------------------------------------------------------
!  Radiation and other stuff

	do j=1,ny
	 do i=1,nx
	   cwp(i,j)=0.
	   cwpl(i,j)=0.
	   cwpm(i,j)=0.
	   cwph(i,j)=0.
	   topind(i,j)=1
           z_inv_ind(i,j)=1
           z_base_ind(i,j)=0
           z_top_ind(i,j)=0
           grad_max(i,j) = 0.
	 end do
	end do
	
	if(CEM) then
	  cwpmax=0.02
	else
	  cwpmax=0.0
	endif

	do k=nzm,1,-1
	 prof1(k)=(radqrlw(k)+radqrsw(k))*factor_xy
	 tmp(1)=rho(k)*adzw(k)*dz
         kc = min(nzm,k+1)
         kb = max(1,k-1)
         tmp(2)=1./(z(kc)-z(kb))
	 do j=1,ny
	  do i=1,nx
            if(z(k).lt.4000.) then ! shallow clouds only
               ! find height of max pot. temp. vert gradient (inversion height) 
               grad = (t(i,j,kc)-t(i,j,kb))*tmp(2)
               if(grad_max(i,j).lt.grad)then
                 grad_max(i,j) = grad
                 z_inv_ind(i,j)=k
               end if
            end if
	    cwp(i,j)=cwp(i,j)+tmp(1)*(qcl(i,j,k)+qci(i,j,k))
            if(pres(k).ge.700.) then
	      cwpl(i,j)=cwpl(i,j)+tmp(1)*(qcl(i,j,k)+qci(i,j,k))
            else if(pres(k).le.400.) then
	      cwph(i,j)=cwph(i,j)+tmp(1)*(qcl(i,j,k)+qci(i,j,k))
            else
	      cwpm(i,j)=cwpm(i,j)+tmp(1)*(qcl(i,j,k)+qci(i,j,k))
	    end if
	    if(cwp(i,j).gt.cwpmax.and.topind(i,j).eq.1)topind(i,j)=k
	  end do
	 end do
	end do

        ncloud = 0
        do k=1,nzm
         do j=1,ny
          do i=1,nx
            if(z_base_ind(i,j).eq.0.and.qcl(i,j,k).gt.0.)then
                z_base_ind(i,j)=k
                ncloud = ncloud+1
            end if
            if(qcl(i,j,k).gt.0.) then
                z_top_ind(i,j)=k
            end if
          end do
         end do
        end do
        if(ncloud.eq.0) then
         coef = 0.
        else
         coef = float(nx*ny)/float(ncloud)
         ncloudy = ncloudy+1
        end if

	do j=1,ny
	 do i=1,nx
	   if(cwp(i,j).gt.cwpmax) s_acld=s_acld+1.
	   if(cwpl(i,j).gt.cwpmax) s_acldl=s_acldl+1.
	   if(cwpm(i,j).gt.cwpmax) s_acldm=s_acldm+1.
	   if(cwph(i,j).gt.cwpmax) s_acldh=s_acldh+1.
	   if(tabs(i,j,topind(i,j)).lt.245.) s_acldcold=s_acldcold+1
           s_sst = s_sst + sstxy(i,j)
           zzz = z(z_inv_ind(i,j))*0.001
           z_inv = z_inv + zzz
           z2_inv = z2_inv + zzz**2
           cwpmean = cwpmean + cwp(i,j)
           cwp2 = cwp2 + cwp(i,j)**2
           if(z_base_ind(i,j).gt.0) then
             z_cbmn = z_cbmn + z(z_base_ind(i,j))*0.001*coef
             z2_cb = z2_cb + (z(z_base_ind(i,j))*0.001)**2*coef
             z_cb = min(z_cb,z(z_base_ind(i,j))*0.001)
           end if
           if(z_top_ind(i,j).gt.0) then
             z_ctmn = z_ctmn + z(z_top_ind(i,j))*0.001*coef
             z2_ct = z2_ct + (z(z_top_ind(i,j))*0.001)**2*coef
             z_ct = max(z_ct,z(z_top_ind(i,j))*0.001)
           end if
	 end do
	end do

        !bloss: compute cumulative cloud fraction going both up and down
        !bloss(2020-11): Two thresholds: cwp_threshold1 = 20 g/m2 (tau~3, assuming effr~10um),
        !                                cwp_threshold2 = 0.2 g/m2 (tau~0.03, assuming effr~10um)
        cwp_threshold1 = 2.e-2
        cldcumup1(:) = 0.
        cldcumdn1(:) = 0.
        cwp_threshold2 = 2.e-4
        cldcumup2(:) = 0.
        cldcumdn2(:) = 0.
        do j = 1,ny
          do i = 1,nx
            ! compute upward cwp, when it hits cwp_threshold1, add one to cldcumup for all levels above
            tmpcwp = 0.
            do k = 1,nzm
              tmpcwp = tmpcwp + rho(k)*adzw(k)*dz*(qcl(i,j,k)+qci(i,j,k))
              if(tmpcwp.gt.cwp_threshold1) then
                cldcumup1(k:nzm) = cldcumup1(k:nzm) + 1.
                EXIT
              end if
            end do
       
            ! compute downward cwp, when it hits cwp_threshold1, add one to cldcum\dn for all levels below
            tmpcwp = 0.
            do k = nzm,1,-1
              tmpcwp = tmpcwp + rho(k)*adzw(k)*dz*(qcl(i,j,k)+qci(i,j,k))
              if(tmpcwp.gt.cwp_threshold1) then
                cldcumdn1(1:k) = cldcumdn1(1:k) + 1.
                EXIT
              end if
            end do

            ! compute upward cwp, when it hits cwp_threshold2, add one to cldcumup for all levels above
            tmpcwp = 0.
            do k = 1,nzm
              tmpcwp = tmpcwp + rho(k)*adzw(k)*dz*(qcl(i,j,k)+qci(i,j,k))
              if(tmpcwp.gt.cwp_threshold2) then
                cldcumup2(k:nzm) = cldcumup2(k:nzm) + 1.
                EXIT
              end if
            end do
       
            ! compute downward cwp, when it hits cwp_threshold2, add one to cldcum\dn for all levels below
            tmpcwp = 0.
            do k = nzm,1,-1
              tmpcwp = tmpcwp + rho(k)*adzw(k)*dz*(qcl(i,j,k)+qci(i,j,k))
              if(tmpcwp.gt.cwp_threshold2) then
                cldcumdn2(1:k) = cldcumdn2(1:k) + 1.
                EXIT
              end if
            end do

          end do
        end do

	call hbuf_put('CLDCUMU1',cldcumup1,factor_xy)
	call hbuf_put('CLDCUMD1',cldcumdn1,factor_xy)
 
	call hbuf_put('CLDCUMU2',cldcumup2,factor_xy)
	call hbuf_put('CLDCUMD2',cldcumdn2,factor_xy)
 
	do k=1,nzm
	 prof2(k)=0.
	 prof3(k)=0.	 
	 n=0
	 if(dolongwave.or.doshortwave) then
	   do j=1,ny
	     do i=1,nx
	       if(cwp(i,j).gt.cwpmax) then
	         n=n+1
	         prof2(k)=prof2(k)+qrad(i,j,k)
	       else
	         prof3(k)=prof3(k)+qrad(i,j,k)
	       endif
	     end do
	   end do 
	 end if
	end do


	call hbuf_put('HLADV',tadv,factor_xy*86400./dtn)
	call hbuf_put('HLDIFF',tdiff,factor_xy*86400./dtn)
	call hbuf_put('HLLAT',tlat+tlatqi,factor_xy*86400./dtn)
	call hbuf_put('HLRAD',prof1,86400.)


        if(dotracers) then
          do ntr=1,ntracers
           call hbuf_avg_put(trim(tracername(ntr)),tracer(:,:,:,ntr), &
                                     dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)
           call hbuf_put(trim(tracername(ntr))//'FLX',trwle(:,ntr),factor_xy)
           call hbuf_put(trim(tracername(ntr))//'FLXS',trwsb(:,ntr),factor_xy)
           call hbuf_put(trim(tracername(ntr))//'ADV',tradv(:,ntr),factor_xy*86400./dtn)
           call hbuf_put(trim(tracername(ntr))//'DIFF',trdiff(:,ntr),factor_xy*86400./dtn)
           call hbuf_put(trim(tracername(ntr))//'PHYS',trphys(:,ntr),factor_xy*86400./dtn)
          end do
        end if


	call hbuf_put('RADLWUP',radlwup,factor_xy)
	call hbuf_put('RADLWDN',radlwdn,factor_xy)
	call hbuf_put('RADSWUP',radswup,factor_xy)
	call hbuf_put('RADSWDN',radswdn,factor_xy)
	call hbuf_put('RADQRLW',radqrlw,factor_xy*86400.) 	
	call hbuf_put('RADQRSW',radqrsw,factor_xy*86400.) 	
	call hbuf_put('RADQR',prof1,86400.) 	
	call hbuf_put('RADQRC',prof2,86400./(n+1.e-5)) 	
	call hbuf_put('RADQRS',prof3,86400./(nx*ny-n+1.e-5))

        if(do_output_clearsky_heating_profiles) then
          call hbuf_put('RADQRCLW',radqrclw,factor_xy*86400.)
          call hbuf_put('RADQRCSW',radqrcsw,factor_xy*86400.)
        end if

        ! Call instrument simulators.  Code is in SRC/SIMULATORS/
        !   As of Feb 2016, this includes the ISCCP, MODIS and MISR
        !   simulators from COSP v1.4.
        call compute_instr_diags()
        ! add output for WTG vertical velocity
        if(dowtg_blossey_etal_JAMES2009) then
          call hbuf_put('WWTG',w_wtg,1.)
          call hbuf_put('TVPR_WTG',tvpr_wtg,1.)
        end if

        ! add output for reference large-scale vertical velocity
        if(dowtg_blossey_etal_JAMES2009) then
          call hbuf_put('WOBSREF',wsub_ref,1.)
        end if

!---------------------------------------------------------
!  Apparent heat/moisture sources/sinks

	tmp(1)=1./dtn

	do k=1,nzm
	 prof2(k)=0.
	 prof3(k)=0.	 
	 n=0
	 do j=1,ny
	  do i=1,nx
	    prof2(k)=prof2(k)+(tabs(i,j,k)-t01(k))*tmp(1)-ttend(k)-prof1(k)
	    prof3(k)=prof3(k)-fac_cond*((qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q01(k))*tmp(1)-qtend(k))
	  end do
	 end do 
	end do

	call hbuf_put('Q1C',prof2,factor_xy*86400.) 	
	call hbuf_put('Q2',prof3,factor_xy*86400.) 	

!===================================================
! UW ADDITIONS

        ! compute storage terms
        ustor(:) = (u0(1:nzm) - ustor(:))/dt/float(nstatis)
        vstor(:) = (v0(1:nzm) - vstor(:))/dt/float(nstatis)
        tstor(:) = (t0(1:nzm) - tstor(:))/dt/float(nstatis)
        qstor(:) = (q0(1:nzm) - qstor(:))/dt/float(nstatis)

	call hbuf_put('USTOR',ustor,86400.)
	call hbuf_put('VSTOR',vstor,86400.)
	call hbuf_put('HLSTOR',tstor,86400.)
	call hbuf_put('QTSTOR',qstor,86400.*1.e3)

        !bloss: extra nudging/vertical large-scale advective tendency outputs
	call hbuf_put('TVTEND',tlsvadv,86400.)
	call hbuf_put('QVTEND',qlsvadv,86400.*1.e3)
	call hbuf_put('THTEND',ttend-tlsvadv,86400.)
	call hbuf_put('QHTEND',qtend-qlsvadv,86400.*1.e3)
	call hbuf_put('TNUDGE',tnudge,86400.)
	call hbuf_put('QNUDGE',qnudge,86400.*1.e3)

	call hbuf_put('ULSADVV',ulsvadv,86400.)
	call hbuf_put('VLSADVV',vlsvadv,86400.)
	call hbuf_put('UNUDGE',unudge,86400.)
	call hbuf_put('VNUDGE',vnudge,86400.)

        !bloss: add outputs for resolved and subgrid momentum tendencies

        ! subgrid tendencies
        udiff(1:nzm) = &
             (rhow(1:nzm)*uwsb(1:nzm) - rhow(2:nz)*uwsb(2:nz)) &
             /rho(1:nzm)/dz/adz(1:nzm)
        vdiff(1:nzm) = &
             (rhow(1:nzm)*vwsb(1:nzm) - rhow(2:nz)*vwsb(2:nz)) &
             /rho(1:nzm)/dz/adz(1:nzm)

        ! resolved tendencies
        uadv(1:nzm) = & 
             (rhow(1:nzm)*uwle(1:nzm) - rhow(2:nz)*uwle(2:nz)) &
             /rho(1:nzm)/dz/adz(1:nzm) - udiff(1:nzm)
        vadv(1:nzm) = & 
             (rhow(1:nzm)*vwle(1:nzm) - rhow(2:nz)*vwle(2:nz)) &
             /rho(1:nzm)/dz/adz(1:nzm) - vdiff(1:nzm)

	call hbuf_put('UADV',uadv,factor_xy*86400.)
	call hbuf_put('VADV',vadv,factor_xy*86400.)
	call hbuf_put('UDIFF',udiff,factor_xy*86400.)
	call hbuf_put('VDIFF',vdiff,factor_xy*86400.)

        ! coriolis accelerations
        call hbuf_put('UTENDCOR',utendcor,factor_xy*86400.)
	call hbuf_put('VTENDCOR',vtendcor,factor_xy*86400.)

        ! momentum budget residual
        call hbuf_put('URESID',ustor-factor_xy*(utendcor+uadv+udiff)-unudge-ulsvadv,1.)
        call hbuf_put('VRESID',vstor-factor_xy*(vtendcor+vadv+vdiff)-vnudge-vlsvadv,1.)

        !bloss: set up storage terms for next time
        ustor(:) = u0(1:nzm)
        vstor(:) = v0(1:nzm)
        tstor(:) = t0(1:nzm)
        qstor(:) = q0(1:nzm)

        if(dosurface) then
          call hbuf_put('SURFVARS',surffluxstats,factor_xy)
        end if

! END UW ADDITIONS
!===================================================

        call t_stopf('statistics')

end
	
