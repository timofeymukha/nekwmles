$(OBJDIR)/spalding.o        :$(NEK_WMLES)/spalding.f;                $(F77) -c $(FL2) -I./ $< -o $@$
$(OBJDIR)/algebraicwm.o     :$(NEK_WMLES)/algebraicwm.f;             $(F77) -c $(FL2) -I./ $< -o $@$
$(OBJDIR)/wmles.o           :$(NEK_WMLES)/wmles.f;                   $(F77) -c $(FL2) -I./ $< -o $@$