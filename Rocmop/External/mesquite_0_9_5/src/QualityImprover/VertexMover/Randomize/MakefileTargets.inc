all_headers: $(RANDSMOOTHHEADERS)
all_objects: $(RANDSMOOTHOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(RANDSMOOTHDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(RANDSMOOTHDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(RANDSMOOTHDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

