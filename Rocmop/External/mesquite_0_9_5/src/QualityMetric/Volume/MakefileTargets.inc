all_headers: $(QMVOLUMEHEADERS)
all_objects: $(QMVOLUMEOBJ)

# puts links to the headers in the include directory    
$(localincludedir)/%.hpp : $(QMVOLUMEDIR)/%.hpp
	${PREFIX} ln -sf $< $@

$(localincludedir)/%.hh : $(QMVOLUMEDIR)/%.hh
	${PREFIX} ln -sf $< $@

$(localobjdir)/%.o : $(QMVOLUMEDIR)/%.cpp
	@echo "Compiling $<"
	${PREFIX} ${CXX} $(CONFIG_CFLAGS) ${CXXFLAGS} -c -o $@ $<

