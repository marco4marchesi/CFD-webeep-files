################################################################################
#
#  UNSTRUCTURED-HYBRID MESH GENERATOR
#
#  Department of Aerospace Engineering
#  Politecnico di Milano. 1995-2008
#  CFD Laboratory, Department of Mechanical Engineering, McGill University, 2009-2015
#  ACE, Department of Mechanical and Aerospace Engineering, University of Strathclyde, 2016-
#
################################################################################

#
SHELL = /bin/sh
UNAME = /bin/uname
OPERSYS = $(shell $(UNAME))
ARCH = $(shell $(UNAME) -m)
#
#----------------------------------------
# COMPILER SETTINGS AND COMMANDS
#----------------------------------------
#
include ../Makefile.config
COMPILE = $(CMP_F90) $(FFLAGS)
LINK = $(LNK_F90) $(FFLAGS)
EXETAG = $(OPERSYS)_$(ARCH)


#----------------------------------------
# LOCAL OBJECTS
#----------------------------------------
OBJ = dynamic_vector.o  banded.o  linsys.o   \
      eigsys_2d.o  grid_types.o  curves.o    \
      list.o  delaunay0.o  eucl_delaunay.o   \
      metric_2d.o  convert.o  hyb_arrays.o   \
      refine_2d.o  back_2d.o  bound_2d.o     \
      check_2d.o  front_2d.o  smooth_bline.o \
      grid_1d.o init_data.o schemi.o  qgraph.o         \
      wakes_gen.o  fork.o  interfer.o        \
      element_topology.o  mesh_structure.o   \
      node_pair_structure.o  nodes.o         \
      np_topology_gen.o            \
      grid_2d.o grid_3d.o import_module.o    \
      export_module.o mesh_transformation.o

#----------------------------------------
# MAIN RULES
#----------------------------------------
all: uhMesh.exe move

uhMesh.exe: uhMesh.o $(OBJ)
	($(LINK) -o uhMesh.$(EXETAG) uhMesh.o $(OBJ))

move:
	mv uhMesh.$(EXETAG) ../bins        

#----------------------------------------
# RULES FOR LOCAL OBJECTS
#----------------------------------------
dynamic_vector.o:
	($(COMPILE) dynamic_vector.f90)
eigsys_2d.o: 
	($(COMPILE) eigsys_2d.f90)
banded.o: 
	($(COMPILE) banded.f90)
linsys.o:
	($(COMPILE) linsys.f90)
grid_types.o:
	($(COMPILE) grid_types.f90)
smooth_bline.o:
	($(COMPILE) smooth_bline.f90)
curves.o: banded.o
	($(COMPILE) curves.f90)
grid_1d.o:
	($(COMPILE) grid_1d.f90)        
init_data.o: curves.o smooth_bline.o grid_1d.o
	($(COMPILE) init_data.f90)
list.o: grid_types.o
	($(COMPILE) list.f90)
delaunay0.o: list.o \
             linsys.o	
	($(COMPILE) delaunay0.f90)
eucl_delaunay.o: delaunay0.o
	($(COMPILE) eucl_delaunay.f90)
wakes_gen.o: list.o
	($(COMPILE) wakes_gen.f90)
fork.o: eucl_delaunay.o \
        grid_types.o    \
        wakes_gen.o     \
        init_data.o  
	($(COMPILE) fork.f90)
interfer.o: grid_types.o \
            eucl_delaunay.o
	($(COMPILE) interfer.f90)
schemi.o: list.o          \
          eucl_delaunay.o \
          init_data.o     \
          grid_types.o    \
          wakes_gen.o     \
          fork.o          \
          interfer.o
	($(COMPILE) schemi.f90)
qgraph.o: grid_types.o \
          schemi.o
	($(COMPILE) qgraph.f90)
convert.o: init_data.o
	($(COMPILE) convert.f90)
metric_2d.o: eigsys_2d.o \
             init_data.o
	($(COMPILE) metric_2d.f90)
hyb_arrays.o: convert.o \
              schemi.o
	($(COMPILE) hyb_arrays.f90)
back_2d.o: hyb_arrays.o \
           metric_2d.o 
	($(COMPILE) back_2d.f90)
bound_2d.o: back_2d.o \
            schemi.o  \
            qgraph.o
	($(COMPILE) bound_2d.f90)
refine_2d.o: bound_2d.o
	($(COMPILE) refine_2d.f90)
check_2d.o:     
	($(COMPILE) check_2d.f90)
front_2d.o: refine_2d.o \
            check_2d.o
	($(COMPILE) front_2d.f90)
element_topology.o: dynamic_vector.o
	($(COMPILE) element_topology.f90)
nodes.o: dynamic_vector.o
	($(COMPILE) nodes.f90)
mesh_structure.o: dynamic_vector.o \
	          nodes.o \
	          element_topology.o
	($(COMPILE) mesh_structure.f90)
node_pair_structure.o: dynamic_vector.o \
	               nodes.o \
	               mesh_structure.o
	($(COMPILE) node_pair_structure.f90)
np_topology_gen.o: dynamic_vector.o \
	           nodes.o \
	           mesh_structure.o \
                   element_topology.o \
                   node_pair_structure.o
	($(COMPILE) np_topology_gen.f90)
grid_2d.o: front_2d.o  \
 	   init_data.o \
           qgraph.o \
           nodes.o
	($(COMPILE) grid_2d.f90)
grid_3d.o: dynamic_vector.o \
           grid_1d.o
	($(COMPILE) grid_3d.f90)
import_module.o: dynamic_vector.o \
	         nodes.o \
 	         mesh_structure.o \
                 node_pair_structure.o \
                 np_topology_gen.o \
                 grid_2d.o \
                 grid_3d.o
	($(COMPILE) import_module.f90)
export_module.o: dynamic_vector.o \
	         nodes.o \
 	         mesh_structure.o \
                 node_pair_structure.o \
                 np_topology_gen.o \
                 grid_1d.o \
                 import_module.o
	($(COMPILE) export_module.f90)
mesh_transformation.o: nodes.o \
                       init_data.o
	($(COMPILE) mesh_transformation.f90)
uhMesh.o: $(OBJ) 
	($(COMPILE) uhMesh.f90)

#-------------------------------------------------------------------------------
#  CLEANING COMMANDS
#-------------------------------------------------------------------------------
.PHONY: cleanLoc cleanAll

clean:
	rm -f ../bins/darwin.* *.o *.mod
