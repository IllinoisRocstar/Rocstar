all_headers: $(QUALIMPROVERHEADERS)
all_objects: $(QUALIMPROVEROBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QUALIMPROVERDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QUALIMPROVERDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QUALIMPROVERDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

