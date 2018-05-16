all_headers: $(FEASNEWTHEADERS)
all_objects: $(FEASNEWTOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(FEASNEWTDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(FEASNEWTDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(FEASNEWTDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

