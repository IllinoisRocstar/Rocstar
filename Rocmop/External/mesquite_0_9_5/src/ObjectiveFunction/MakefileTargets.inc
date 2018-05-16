all_headers: $(OBJFUNCHEADERS)
all_objects: $(OBJFUNCOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(OBJFUNCDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(OBJFUNCDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(OBJFUNCDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

