all_headers: $(QMUNTANGLEHEADERS)
all_objects: $(QMUNTANGLEOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMUNTANGLEDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QMUNTANGLEDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QMUNTANGLEDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

