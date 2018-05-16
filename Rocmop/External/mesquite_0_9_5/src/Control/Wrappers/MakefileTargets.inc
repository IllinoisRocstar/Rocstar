all_headers: $(WRAPPERSHEADERS)
all_objects: $(WRAPPERSOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(WRAPPERSDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(WRAPPERSDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(WRAPPERSDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

