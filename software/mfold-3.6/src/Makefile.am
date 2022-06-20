AM_CPPFLAGS = -DPKGDATADIR=\"$(pkgdatadir)/\"
AM_FFLAGS = -C -O3 -Wall -fno-automatic -ffixed-line-length-none 

EXTRA_DIST = pre-efiles.f pre-main.f pre-newtemp.f pre-efn.f \
pre-misc-quik.f

CLEANFILES = efiles.f efn.f misc-quik.f newtemp.f main.f
DISTCLEANFILES = efiles.f efn.f misc-quik.f newtemp.f main.f

bin_PROGRAMS = distance ct2bp ct_compare efn newtemp scorer auxgen \
sav2plot sav2p_num h_num ss_count add_dHdSTm nafold quikfold \
add_dHdSTm2 boxplot_ng ct_boxplot_ng overlay_boxplot_ng sir_graph mfold_datdir

mfold_datdir_SOURCES = mfold_datdir.c

distance_SOURCES = distance.f bpin.f bpsrch.f
distance_FFLAGS = $(AM_FFLAGS)

ct_compare_SOURCES = ct_compare.f
ct_compare_FFLAGS = $(AM_FFLAGS)

ct2bp_SOURCES = ct2bp.f
ct2bp_FFLAGS = $(AM_FFLAGS)

efn_SOURCES = efn.f efn.inc
efn_FFLAGS = $(AM_FFLAGS)

h_num_SOURCES = h-num.f
h_num_FFLAGS = $(AM_FFLAGS)

newtemp_SOURCES = newtemp.f ion.f dc.f 
newtemp_FFLAGS = $(AM_FFLAGS) 

scorer_SOURCES = scorer.f
scorer_FFLAGS = $(AM_FFLAGS)

auxgen_SOURCES = auxgen.f
auxgen_FFLAGS = $(AM_FFLAGS)

sav2plot_SOURCES = sav2plot.f maxn.inc
sav2plot_FFLAGS = $(AM_FFLAGS)

sav2p_num_SOURCES = sav2p-num.f maxn.inc
sav2p_num_FFLAGS = $(AM_FFLAGS)

ss_count_SOURCES = ss-count.f
ss_count_FFLAGS = $(AM_FFLAGS)

add_dHdSTm_SOURCES = add-dHdSTm.f
add_dHdSTm_FFLAGS = $(AM_FFLAGS)

add_dHdSTm2_SOURCES = add-dHdSTm2.f
add_dHdSTm2_FFLAGS = $(AM_FFLAGS)

boxplot_ng_SOURCES = boxplot_ng.c boxplot_setcolor.h boxplot_input.h \
boxplot_ps.h boxplot_img.h strings_img.h util.h

ct_boxplot_ng_SOURCES = ct_boxplot_ng.c ct_boxplot_general.h \
ct_boxplot_img.h ct_boxplot_ps.h ct_boxplot_read_ct.h \
ct_boxplot_setcolor.h ct_boxplot.h header_ps.h strings_img.h util.h

overlay_boxplot_ng_SOURCES = overlay_boxplot_ng.c overlay_boxplot.h \
overlay_boxplot_setcolor.h overlay_boxplot_read_ct.h strings_img.h \
overlay_boxplot_ps.h overlay_boxplot_img.h header_ps.h util.h

sir_graph_SOURCES = sir_graph.c sir_graph_color.h sir_graph_color_ann.h \
sir_graph_ann.h sir_graph_ps.h sir_graph_img.h sir_graph_ss.h sir_graph_ct.h \
sir_graph_alloc.h strings_img.h util.h

BITS = bits-gnu.f

INCL = rna.inc

NASRC = main.f rna.f misc.f efiles.f sort.f dc.f $(BITS) formid.f multid.f

QUIK = quik.f rna-quik.f misc-quik.f dc.f multid.f

nafold_SOURCES = $(INCL) maxn.inc $(NASRC)
nafold_FFLAGS = $(AM_FFLAGS)

quikfold_SOURCES = maxn1.inc quik.inc $(QUIK)
quikfold_FFLAGS = $(AM_FFLAGS)

newtemp.f: pre-newtemp.f
	@@AWK@ -v ddir=$(pkgdatadir) '{sub("PKGDATADIR",ddir); print $0}' pre-newtemp.f > newtemp.f

efiles.f: pre-efiles.f
	@@AWK@ -v ddir=$(pkgdatadir) '{sub("PKGDATADIR",ddir); print $0}' pre-efiles.f > efiles.f

efn.f: pre-efn.f
	@@AWK@ -v ddir=$(pkgdatadir) '{sub("PKGDATADIR",ddir); print $0}' pre-efn.f > efn.f

main.f: pre-main.f
	@@AWK@ -v pack=$(PACKAGE) -v vers=$(VERSION) \
	'{sub("PACKAGE",pack); sub("VERSION",vers); print $0}' \
	pre-main.f > main.f

misc-quik.f: pre-misc-quik.f
	@@AWK@ -v ddir=$(pkgdatadir) '{sub("PKGDATADIR",ddir); print $0}' pre-misc-quik.f > misc-quik.f

