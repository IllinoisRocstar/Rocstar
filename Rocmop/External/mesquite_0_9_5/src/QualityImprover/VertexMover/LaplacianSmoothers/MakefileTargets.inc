all_headers: $(LAPLSMOOTHHEADERS)
all_objects: $(LAPLSMOOTHOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(LAPLSMOOTHDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(LAPLSMOOTHDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(LAPLSMOOTHDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

