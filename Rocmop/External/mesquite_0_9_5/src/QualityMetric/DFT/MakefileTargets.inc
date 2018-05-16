all_headers: $(DFTHEADERS)
all_objects: $(DFTOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(DFTDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(DFTDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(DFTDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

