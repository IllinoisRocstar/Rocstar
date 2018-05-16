all_headers: $(STEEPDESCHEADERS)
all_objects: $(STEEPDESCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(STEEPDESCDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(STEEPDESCDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(STEEPDESCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

