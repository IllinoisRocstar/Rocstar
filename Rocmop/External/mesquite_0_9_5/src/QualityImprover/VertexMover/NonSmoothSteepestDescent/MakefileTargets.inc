all_headers: $(NONSMOOTHSTEEPDESCHEADERS)
all_objects: $(NONSMOOTHSTEEPDESCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(NONSMOOTHSTEEPDESCDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(NONSMOOTHSTEEPDESCDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(NONSMOOTHSTEEPDESCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -DMSQ_DBG0 -c -o $@ $<

