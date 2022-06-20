bin_SCRIPTS = h-num mfold mfold_quik reformat-seq.sh filter-sort myps2img.bash

EXTRA_DIST = myps2img.bash

CLEANFILES = mfold mfold_quik h-num reformat-seq filter-sort 
DISTCLEANFILES = $(CLEANFILES)

install-exec-hook :
	cd $(DESTDIR)$(bindir) && $(LN_S) -f myps2img.bash myps2pdf.bash
	cd $(DESTDIR)$(bindir) && $(LN_S) -f myps2img.bash myps2jpg.bash
	cd $(DESTDIR)$(bindir) && $(LN_S) -f myps2img.bash myps2png.bash

