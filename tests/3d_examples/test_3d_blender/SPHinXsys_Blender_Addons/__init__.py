# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

bl_info = {
    "name" : "SPHinXsys",
    "author" : "Chi Zhang & Xiangyu Hu",
    "description" : "SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle Hydrodynamics for industrial compleX systems. It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization.",
    "blender" : (2, 80, 0),
    "version" : (0, 0, 1),
    "location" : "",
    "warning" : "",
    "category" : "Generic"
}

import bpy
from .props import SPHinXsysGlobalProps, SPHinXsysProps, registerSPHinXsysProperties
from .saveMesh import Mesh_OT_SaveSelectdMeshes2File, Mesh_OT_SaveAllMeshes2File, Mesh_OT_ClearSTLFiles
from .panel import SPHinXsys_PT_Panel
from .SPHinXsys import SPHinXsys_OT_runCase

classes = [SPHinXsysGlobalProps, SPHinXsysProps, Mesh_OT_ClearSTLFiles, Mesh_OT_SaveSelectdMeshes2File, Mesh_OT_SaveAllMeshes2File, SPHinXsys_PT_Panel, SPHinXsys_OT_runCase]
#register, unregister = bpy.utils.register_classes_factory(classes)

def register():
    for item in classes:
        bpy.utils.register_class(item)

    registerSPHinXsysProperties()

def unregister():
    for item in classes:
        bpy.utils.unregister_class(item)



