all_headers: $(MISCHEADERS)
all_objects: $(MISCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(MISCDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(MISCDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(MISCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

