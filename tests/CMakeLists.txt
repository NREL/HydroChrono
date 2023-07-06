find_package(Python COMPONENTS Interpreter REQUIRED)

add_executable(h5fileinfo_t01 h5fileinfo_t01.cpp)
target_link_libraries(h5fileinfo_t01 HydroChrono)

add_executable(chloadaddedmass_t01 chloadaddedmass_t01.cpp)
target_link_libraries(chloadaddedmass_t01 HydroChrono)

add_executable(chrono_error_t01 chrono_error_t01.cpp)
target_link_libraries(chrono_error_t01 HydroChrono)

# ============
# TESTS
# ============
#ex: to runn all test with label  demos AND label small 
# > ctest -L demos -L small

if(TARGET h5fileinfo_t01)
        add_test (
                NAME h5fileinfo_01
                COMMAND $<TARGET_FILE:h5fileinfo_t01> ${HYDROCHRONO_DATA_DIR}
        )
        set_tests_properties(
            h5fileinfo_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET h5fileinfo_t01)

if(TARGET chloadaddedmass_t01)
        add_test (
                NAME chloadaddedmass_01
                COMMAND $<TARGET_FILE:chloadaddedmass_t01> ${HYDROCHRONO_DATA_DIR} 
        )
        set_tests_properties(
                chloadaddedmass_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET chloadaddedmass_t01)


if(TARGET chrono_error_t01)
        add_test (
                NAME chrono_error_01
                COMMAND $<TARGET_FILE:chrono_error_t01> ${HYDROCHRONO_DATA_DIR}
        )
        set_tests_properties(
                chrono_error_01
                PROPERTIES LABELS "examples;small;core"
        )
endif(TARGET chrono_error_t01)

# DEMO SPHERE


# TODO : ensure previous result test file is removed  before comparing results
# The fixture dependency is not sufficient. It forces relaunche of the case but 
# in case of failure the previous file is still there for comparison ! 
if(TARGET demo_sphere_decay)


        add_test (
                NAME demo_sphere_decay_01
                COMMAND $<TARGET_FILE:demo_sphere_decay> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        
        set_tests_properties(
                demo_sphere_decay_01
                PROPERTIES 
                LABELS "demos;small;core"
                FIXTURES_SETUP  demo_sphere_decay_01_file
        )

        set(FILE_RST ${CMAKE_CURRENT_BINARY_DIR}/results/sphere_decay.txt)
        set(FILE_REF ${HYDROCHRONO_DATA_DIR}/sphere/postprocessing/sphere_decay_comparison.txt)

        add_test(NAME demo_sphere_decay_ref
                COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare.py ${FILE_REF} ${FILE_RST}
        )

        set_tests_properties(
                demo_sphere_decay_ref
                PROPERTIES 
                LABELS "demos;ref"
                DEPENDS demo_sphere_decay_01
                REQUIRED_FILES ${FILE_RST} 
                FIXTURES_REQUIRED demo_sphere_decay_01_file
        )

endif()

# TODO Add test  comparisons
if(TARGET demo_sphere_reg_waves)
        add_test (
                NAME demo_sphere_reg_waves_01
                COMMAND $<TARGET_FILE:demo_sphere_reg_waves> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_sphere_reg_waves_01
                PROPERTIES LABELS "demos;long;core"
        )
endif()


if(TARGET demo_sphere_irreg_waves)
        add_test (
                NAME demo_sphere_irreg_waves_01
                COMMAND $<TARGET_FILE:demo_sphere_irreg_waves> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_sphere_irreg_waves_01
                PROPERTIES LABELS "medium;core"
        )
endif()


# DEMO RM3
if(TARGET demo_rm3_decay)
        add_test (
                NAME demo_rm3_decay_01
                COMMAND $<TARGET_FILE:demo_rm3_decay> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_rm3_decay_01
                PROPERTIES LABELS "demos;medium;core"
        )
endif() 

if(TARGET demo_rm3_reg_waves)
        add_test (
                NAME demo_rm3_reg_waves_01
                COMMAND $<TARGET_FILE:demo_rm3_reg_waves> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_rm3_reg_waves_01
                PROPERTIES LABELS "demos;long;core"
        )
endif() 

# DEMO OSWEC
if(TARGET demo_oswec_decay)
        add_test (
                NAME demo_oswec_decay_01
                COMMAND $<TARGET_FILE:demo_oswec_decay> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_oswec_decay_01
                PROPERTIES LABELS "demos;medium;core"
        )
endif() 

if(TARGET demo_oswec_reg_waves)
        add_test (
                NAME demo_oswec_reg_waves_01
                COMMAND $<TARGET_FILE:demo_oswec_reg_waves> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_oswec_reg_waves_01
                PROPERTIES LABELS "demos;long;core"
        )
endif() 

# DEMO DeepCWind
if(TARGET demo_DeepCWind_decay)
        add_test (
                NAME demo_DeepCWind_decay_01
                COMMAND $<TARGET_FILE:demo_DeepCWind_decay> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_DeepCWind_decay_01
                PROPERTIES LABELS "demos;medium;core"
        )
endif() 

# DEMO F3OF
if(TARGET demo_f3of_DT1)
        add_test (
                NAME demo_F3OF_DT1_01
                COMMAND $<TARGET_FILE:demo_f3of_DT1> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_F3OF_DT1_01
                PROPERTIES LABELS "demos;long;core"
        )
endif() 

if(TARGET demo_f3of_DT2)
        add_test (
                NAME demo_F3OF_DT2_01
                COMMAND $<TARGET_FILE:demo_f3of_DT2> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_F3OF_DT2_01
                PROPERTIES LABELS "demos;long;core"
        )
endif() 

if(TARGET demo_f3of_DT3)
        add_test (
                NAME demo_F3OF_DT3_01
                COMMAND $<TARGET_FILE:demo_f3of_DT3> ${HYDROCHRONO_DATA_DIR} --nogui
        )
        set_tests_properties(
                demo_F3OF_DT3_01
                PROPERTIES LABELS "demos;long;core"
        )
endif() 



