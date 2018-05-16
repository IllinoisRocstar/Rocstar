all_headers: $(QUALMETRICHEADERS)
all_objects: $(QUALMETRICOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QUALMETRICDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QUALMETRICDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QUALMETRICDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<
