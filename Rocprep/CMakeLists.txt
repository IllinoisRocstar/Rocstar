cmake_minimum_required(VERSION 3.0)

set(PREP_FILES
    ModuleProcessor.pm
    Properties.pm
    RocprepProperties.pm
    RocfloProcessor.pm
    RocfluProcessor.pm
    RocfaceProcessor.pm
    RocfracProcessor.pm
    RocsolidProcessor.pm
    RocstarProcessor.pm
    RunLogControl.pm
    Rocprep.pm
    rocprep
)

file(COPY ${PREP_FILES}
    DESTINATION ${CMAKE_BINARY_DIR}/bin
    FILE_PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)

install(FILES ${PREP_FILES}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
)
