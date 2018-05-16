all_headers: $(CONJGRADHEADERS)
all_objects: $(CONJGRADOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(CONJGRADDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(CONJGRADDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(CONJGRADDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

