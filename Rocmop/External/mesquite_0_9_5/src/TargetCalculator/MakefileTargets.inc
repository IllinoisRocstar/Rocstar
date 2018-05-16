all_headers: $(TARGCALCHEADERS)
all_objects: $(TARGCALCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(TARGCALCDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(TARGCALCDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(TARGCALCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

