all_headers: $(TOPOMODIFHEADERS)
all_objects: $(TOPOMODIFOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(TOPOMODIFDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(TOPOMODIFDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(TOPOMODIFDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

