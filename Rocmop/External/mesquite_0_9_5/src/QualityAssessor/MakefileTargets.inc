all_headers: $(QUALASSESSHEADERS)
all_objects: $(QUALASSESSOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QUALASSESSDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QUALASSESSDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QUALASSESSDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

