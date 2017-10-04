noinst_LIBRARIES = libAtoms.a
libAtoms_a_SOURCES = Atoms.c Atoms.h CHECK.c ChainofStates.c Config.c Dacapo.c Gaussian.c Gr.c LeastSquareStrain.c Mendeleyev.c Motion.c NIST.h Neighborlist.c VASP.c Voronoi.c


AM_CFLAGS = -I$(top_srcdir)/VecMat3 -I$(top_srcdir)/VecMat -I$(top_srcdir)/IO -I$(top_srcdir)/Scalar



generic_LDADD = libAtoms.a $(top_srcdir)/VecMat3/libVecMat3.a $(top_srcdir)/VecMat/libVecMat.a $(top_srcdir)/Scalar/libScalar.a $(top_srcdir)/IO/libIO.a $(top_srcdir)/Timer/libTimer.a -lgsl -lgslcblas -lreadline -lm  -lncurses -lpthread -llapack -lblas

bin_PROGRAMS = mul iniT punchcylinder vcut trans scut voronoirize label_immobile c2p p2c cfg2xyz xtal vasp2cfg vasp2out cfg2vasp chemical_randomize perturb_s perturb_x fch2cfg chg2xsf cfg2dacapo annotate_atomic_strain affine_deform linear_path period_reduce infectH concat div mobilize_all transmute chg_planar_avg


mul_SOURCES = Config.c Neighborlist.c
mul_CFLAGS = -D_mul $(AM_CFLAGS)
mul_LDADD = $(generic_LDADD)


iniT_SOURCES = Motion.c
iniT_CFLAGS = -D_iniT $(AM_CFLAGS)
iniT_LDADD = $(generic_LDADD)


punchcylinder_SOURCES = Config.c
punchcylinder_CFLAGS = -D_punchcylinder $(AM_CFLAGS)
punchcylinder_LDADD = $(generic_LDADD)


vcut_SOURCES = Config.c
vcut_CFLAGS = -D_vcut $(AM_CFLAGS)
vcut_LDADD = $(generic_LDADD)


trans_SOURCES = Config.c Motion.c
trans_CFLAGS = -D_trans $(AM_CFLAGS)
trans_LDADD = $(generic_LDADD)


scut_SOURCES = Config.c
scut_CFLAGS = -D_scut $(AM_CFLAGS)
scut_LDADD = $(generic_LDADD)


voronoirize_SOURCES = Voronoi.c
voronoirize_CFLAGS = -D_voronoirize $(AM_CFLAGS)
voronoirize_LDADD = $(generic_LDADD)


label_immobile_SOURCES = Motion.c
label_immobile_CFLAGS = -D_label_immobile $(AM_CFLAGS)
label_immobile_LDADD = $(generic_LDADD)


c2p_SOURCES = Config.c
c2p_CFLAGS = -D_c2p $(AM_CFLAGS)
c2p_LDADD = $(generic_LDADD)


p2c_SOURCES = Config.c
p2c_CFLAGS = -D_p2c $(AM_CFLAGS)
p2c_LDADD = $(generic_LDADD)


cfg2xyz_SOURCES = Config.c
cfg2xyz_CFLAGS = -D_cfg2xyz $(AM_CFLAGS)
cfg2xyz_LDADD = $(generic_LDADD)


xtal_SOURCES = Config.c
xtal_CFLAGS = -D_xtal $(AM_CFLAGS)
xtal_LDADD = $(generic_LDADD)


vasp2cfg_SOURCES = VASP.c
vasp2cfg_CFLAGS = -D_vasp2cfg $(AM_CFLAGS)
vasp2cfg_LDADD = $(generic_LDADD)


vasp2out_SOURCES = VASP.c
vasp2out_CFLAGS = -D_vasp2out $(AM_CFLAGS)
vasp2out_LDADD = $(generic_LDADD)


cfg2vasp_SOURCES = VASP.c
cfg2vasp_CFLAGS = -D_cfg2vasp $(AM_CFLAGS)
cfg2vasp_LDADD = $(generic_LDADD)


chemical_randomize_SOURCES = Config.c
chemical_randomize_CFLAGS = -D_chemical_randomize $(AM_CFLAGS)
chemical_randomize_LDADD = $(generic_LDADD)


perturb_s_SOURCES = Config.c
perturb_s_CFLAGS = -D_perturb_s $(AM_CFLAGS)
perturb_s_LDADD = $(generic_LDADD)


perturb_x_SOURCES = Config.c
perturb_x_CFLAGS = -D_perturb_x $(AM_CFLAGS)
perturb_x_LDADD = $(generic_LDADD)


fch2cfg_SOURCES = Gaussian.c
fch2cfg_CFLAGS = -D_fch2cfg $(AM_CFLAGS)
fch2cfg_LDADD = $(generic_LDADD)


chg2xsf_SOURCES = VASP.c
chg2xsf_CFLAGS = -D_chg2xsf $(AM_CFLAGS)
chg2xsf_LDADD = $(generic_LDADD)


cfg2dacapo_SOURCES = Dacapo.c
cfg2dacapo_CFLAGS = -D_cfg2dacapo $(AM_CFLAGS)
cfg2dacapo_LDADD = $(generic_LDADD)


annotate_atomic_strain_SOURCES = LeastSquareStrain.c
annotate_atomic_strain_CFLAGS = -D_annotate_atomic_strain $(AM_CFLAGS)
annotate_atomic_strain_LDADD = $(generic_LDADD)


affine_deform_SOURCES = Config.c
affine_deform_CFLAGS = -D_affine_deform $(AM_CFLAGS)
affine_deform_LDADD = $(generic_LDADD)


linear_path_SOURCES = ChainofStates.c
linear_path_CFLAGS = -D_linear_path $(AM_CFLAGS)
linear_path_LDADD = $(generic_LDADD)


period_reduce_SOURCES = Config.c
period_reduce_CFLAGS = -D_period_reduce $(AM_CFLAGS)
period_reduce_LDADD = $(generic_LDADD)


infectH_SOURCES = Config.c
infectH_CFLAGS = -D_infectH $(AM_CFLAGS)
infectH_LDADD = $(generic_LDADD)


concat_SOURCES = Config.c
concat_CFLAGS = -D_concat $(AM_CFLAGS)
concat_LDADD = $(generic_LDADD)


div_SOURCES = Config.c
div_CFLAGS = -D_div $(AM_CFLAGS)
div_LDADD = $(generic_LDADD)


mobilize_all_SOURCES = Motion.c
mobilize_all_CFLAGS = -D_mobilize_all $(AM_CFLAGS)
mobilize_all_LDADD = $(generic_LDADD)


transmute_SOURCES = Motion.c
transmute_CFLAGS = -D_transmute $(AM_CFLAGS)
transmute_LDADD = $(generic_LDADD)


chg_planar_avg_SOURCES = VASP.c
chg_planar_avg_CFLAGS = -D_chg_planar_avg $(AM_CFLAGS)
chg_planar_avg_LDADD = $(generic_LDADD)


