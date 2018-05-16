all_headers: $(QMSHAPEHEADERS)
all_objects: $(QMSHAPEOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMSHAPEDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QMSHAPEDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QMSHAPEDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

