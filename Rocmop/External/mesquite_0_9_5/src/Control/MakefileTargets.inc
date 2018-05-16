all_headers: $(CONTROLHEADERS)
all_objects: $(CONTROLOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(CONTROLDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(CONTROLDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(CONTROLDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

