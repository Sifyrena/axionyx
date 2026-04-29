#ifdef AXIONYX
#include <Nyx.H>

using namespace amrex;

#ifdef FDM
void
Nyx::compute_axionyx_quantities (Real& mass, Real& epot, Real& ekinrho, Real& ekinv, Real& ekin, Real& etot, Real& angmom_x,
				 Real& angmom_y, Real& angmom_z, Real& max_dens)
{
  int             finest_level = parent->finestLevel();
  Real            time         = state[Axion_Type].curTime();

  for (int lev = 0; lev <= finest_level; lev++)
    {
      Nyx& nyx_lev = get_level(lev);
      mass += nyx_lev.vol_weight_sum("AxDens", time,true);
      epot += nyx_lev.vol_weight_sum("fdm_epot", time,true);
      ekinrho += nyx_lev.vol_weight_sum("fdm_ekinrho", time,true);
      ekinv += nyx_lev.vol_weight_sum("fdm_ekinv", time,true);
      ekin += nyx_lev.vol_weight_sum("fdm_ekin", time,true);
      angmom_x += nyx_lev.vol_weight_sum("fdm_angmomx", time,true);
      angmom_y += nyx_lev.vol_weight_sum("fdm_angmomy", time,true);
      angmom_z += nyx_lev.vol_weight_sum("fdm_angmomz", time,true);
    }
  etot = epot + ekin;
  for (int lev = 0; lev <= parent->finestLevel(); lev++)
    max_dens = std::max(max_dens, get_level(lev).get_new_data(Axion_Type).max(Nyx::AxDens));
}

#ifndef NO_HYDRO // YES_HYDRO

void
Nyx::compute_baryon_quantities (Real& Bmass, Real& Bepot, Real& Bekin, Real& Beint, Real& Betot)
{
  int             finest_level = parent->finestLevel();
  Real            time         = state[State_Type].curTime();

  for (int lev = 0; lev <= finest_level; lev++)
    {
      Nyx& nyx_lev = get_level(lev);
      Bmass += nyx_lev.vol_weight_sum("density", time,true);
      Bepot += nyx_lev.vol_weight_sum("poteng", time,true);
      Bekin += nyx_lev.vol_weight_sum("kineng", time,true);
      Beint += nyx_lev.vol_weight_sum("eint_e", time,true);
    }
  Betot = Bepot + Bekin + Beint;
}

#endif

void
Nyx::findCentralFDMHaloPosition (){
#ifdef GRAVITY
  IntVect maxdist(16);
  for(int lev=parent->finestLevel();lev>=0;lev--){
    MultiFab& phi = get_level(lev).get_new_data(PhiGrav_Type);
    const Geometry& geomlev = parent->Geom(lev);
    const Real* dx = geomlev.CellSize();
    IntVect fhp;
    AMREX_D_TERM(fhp[0]=floor((fdm_halo_pos[0] - geomlev.ProbLo(0))/dx[0]);,
		 fhp[1]=floor((fdm_halo_pos[1] - geomlev.ProbLo(1))/dx[1]);,
		 fhp[2]=floor((fdm_halo_pos[2] - geomlev.ProbLo(2))/dx[2]););
    Box b = Box(fhp-maxdist, fhp+maxdist);
    Real phimax = phi.max(b,0);
    if(phimax>std::numeric_limits<Real>::lowest()){
      //fhp = phi.maxIndex(b,0);
      fhp = phi.maxIndex(0);
      fdm_halo_pos[0] = (fhp[0]+0.5)*dx[0]+geomlev.ProbLo(0);
      fdm_halo_pos[1] = (fhp[1]+0.5)*dx[1]+geomlev.ProbLo(1);
      fdm_halo_pos[2] = (fhp[2]+0.5)*dx[2]+geomlev.ProbLo(2);
      break;
    }}
#else
  IntVect maxdist(16);
  for(int lev=parent->finestLevel();lev>=0;lev--){
    MultiFab& Ax_new = get_level(lev).get_new_data(Axion_Type);
    const Geometry& geomlev = parent->Geom(lev);
    const Real* dx = geomlev.CellSize();
    IntVect fhp;
    AMREX_D_TERM(fhp[0]=floor((fdm_halo_pos[0] - geomlev.ProbLo(0))/dx[0]);,
		 fhp[1]=floor((fdm_halo_pos[1] - geomlev.ProbLo(1))/dx[1]);,
		 fhp[2]=floor((fdm_halo_pos[2] - geomlev.ProbLo(2))/dx[2]););
    Box b = Box(fhp-maxdist, fhp+maxdist);
    Real densmax = Ax_new.max(b,Nyx::AxDens);
    if(densmax>std::numeric_limits<Real>::lowest()){
      //fhp = Ax_new.maxIndex(b,Nyx::AxDens);
      fhp = Ax_new.maxIndex(Nyx::AxDens);
      fdm_halo_pos[0] = (fhp[0]+0.5)*dx[0]+geomlev.ProbLo(0);
      fdm_halo_pos[1] = (fhp[1]+0.5)*dx[1]+geomlev.ProbLo(1);
      fdm_halo_pos[2] = (fhp[2]+0.5)*dx[2]+geomlev.ProbLo(2);
      break;
    }}
#endif
  std::string amr_prefix = "amr";
  ParmParse ppamr(amr_prefix);
  Vector<std::string> refinement_indicators;
  ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
  for (int i=0; i<refinement_indicators.size(); ++i){
    std::string ref_prefix = amr_prefix + "." + refinement_indicators[i];
    ParmParse ppr(ref_prefix);
    RealBox realbox;
    if (ppr.countval("fdm_halo_box")){
      std::vector<Real> halo_box(BL_SPACEDIM), box_lo(BL_SPACEDIM), box_hi(BL_SPACEDIM);
      ppr.getarr("fdm_halo_box",halo_box,0,halo_box.size());
      for (int d=0; d < AMREX_SPACEDIM; ++d){
	box_lo[d] = Nyx::fdm_halo_pos[d]-halo_box[d]/2.0;
	box_hi[d] = Nyx::fdm_halo_pos[d]+halo_box[d]/2.0;
      }
      
      realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
      AMRErrorTagInfo info;
        if (realbox.ok()) {
          info.SetRealBox(realbox);
        }
    }
  }
}
#endif

#ifdef AXCOMPLEX
void
Nyx::compute_axionyx_quantities (Real& stringdens, Vector<Vector<Real> >& test_val)
{
  if( (Nyx::nstep_spectrum!=-1) && (parent->levelSteps(0)%Nyx::nstep_spectrum == 0) )
    prepare_and_compute_powerspectrum(state[Axion_Type].curTime());
  stringdens = get_level(parent->finestLevel()).vol_weight_sum("axion_string", state[Axion_Type].curTime(), false)
    *pow(state[Axion_Type].curTime()/get_level(parent->finestLevel()).geom.CellSizeArray()[0],2)/6.0;
  if (stringdens == 0.0 && parent->finestLevel()>0)
    stringdens = get_level(parent->finestLevel()-1).vol_weight_sum("axion_string", state[Axion_Type].curTime(), false)
      *pow(state[Axion_Type].curTime()/get_level(parent->finestLevel()-1).geom.CellSizeArray()[0],2)/6.0;

  MultiFab& Ax_new = get_level(0).get_new_data(Axion_Type);
  test_val.resize(Nyx::test_posx.size());
  for(int i=0; i<Nyx::test_posx.size(); i++){
    test_val[i].resize(Ax_new.nComp());
    for (MFIter mfi(Ax_new,false); mfi.isValid(); ++mfi){
      const Box& bx = mfi.validbox();
      const Dim3 lo = amrex::lbound(bx);
      const Dim3 hi = amrex::ubound(bx);
      if(lo.x<=Nyx::test_posx[i] && lo.y<=Nyx::test_posy[i] && lo.z<=Nyx::test_posz[i] && Nyx::test_posx[i]<=hi.x && Nyx::test_posy[i]<=hi.y && Nyx::test_posz[i]<=hi.z){
	for(int n=0; n<Ax_new.nComp(); n++)
	  test_val[i][n] = Ax_new.array(mfi)(Nyx::test_posx[i],Nyx::test_posy[i],Nyx::test_posz[i],n);
	break;
      }}
    ParallelDescriptor::ReduceRealSum(test_val[i].dataPtr(), test_val[i].size(), ParallelDescriptor::IOProcessorNumber());
  }
}
#endif
#endif
