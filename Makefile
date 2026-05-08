MODULE_TOPDIR = ../..

PGM = i.hyper.poppy

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make

default: script html $(TEST_DST)

$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html
