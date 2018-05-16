all_headers: $(QMSMOOTHNESSHEADERS)
all_objects: $(QMSMOOTHNESSOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMSMOOTHNESSDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QMSMOOTHNESSDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QMSMOOTHNESSDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

