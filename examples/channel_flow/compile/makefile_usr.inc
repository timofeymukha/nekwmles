$(OBJDIR)/frame.o           :$(FRAME_SRC)/frame.f;                   $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog_block.o   :$(MONITOR_SRC)/mntrlog_block.f;         $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog.o         :$(MONITOR_SRC)/mntrlog.f;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr_block.o   :$(MONITOR_SRC)/mntrtmr_block.f;         $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr.o         :$(MONITOR_SRC)/mntrtmr.f;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm_block.o      :$(RUNPARAM_SRC)/rprm_block.f;           $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm.o            :$(RUNPARAM_SRC)/rprm.f;                 $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools_block.o  :$(IO_TOOLS_SRC)/io_tools_block.f;       $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools.o        :$(IO_TOOLS_SRC)/io_tools.f;             $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/chkpoint.o        :$(CHKPT_SRC)/chkpoint.f;                $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/chkpt_mstp.o      :$(CHKPTMS_SRC)/chkpt_mstp.f;            $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/map2D.o           :${MAP2D_SRC}/map2D.f;                   $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/stat.o            :./stat.f ./STATD;             $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/stat_IO.o         :$(STAT_SRC)/stat_IO.f ./STATD;          $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/math_tools.o      :$(MATH_SRC)/math_tools.f;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/spalding.o        :$(NEK_WMLES)/spalding.f;                $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/newton.o          :$(NEK_WMLES)/newton.f;                  $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/algebraicwm.o     :$(NEK_WMLES)/algebraicwm.f;             $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/wmles.o           :$(NEK_WMLES)/wmles.f;                   $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/vreman.o          :sgs_models/vreman.f SGS;                $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/utils.o           :sgs_models/utils.f SGS;                 $(F77) -c $(FL2) -I./ $< -o $@
