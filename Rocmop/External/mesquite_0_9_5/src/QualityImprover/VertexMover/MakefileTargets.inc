all_headers: $(VERTEXMOVERHEADERS)
all_objects: $(VERTEXMOVEROBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(VERTEXMOVERDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(VERTEXMOVERDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(VERTEXMOVERDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

