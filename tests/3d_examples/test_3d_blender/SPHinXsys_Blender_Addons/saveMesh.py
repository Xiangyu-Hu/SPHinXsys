import bpy
import os
import sys
projectPath = os.path.abspath('../SPHinXsys_Blender_Addons')
sys.path.append(projectPath)
import props

class Mesh_OT_ClearSTLFiles(bpy.types.Operator):
    bl_idname = "mesh.clearstlfiles"
    bl_label = "Clear STL Files"
    bl_options = {'REGISTER'}

    def execute(self, context):
        props.clearAllBlocksName()
        # props.printBlocksName()
        basedir = os.path.dirname(bpy.data.filepath)
        if not basedir:
            raise Exception("Blend file is not saved")
        basedir += '/SPHinXsys_Blender_Addons/Mesh/'
        fluidsDir = basedir + 'Fluids'
        rigidsDir = basedir + 'Rigids/'
        elasticsDir = basedir + 'Elastics/'

        del_list = os.listdir(fluidsDir)
        for f in del_list:
            file_path = os.path.join(fluidsDir, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
        del_list = os.listdir(rigidsDir)
        for f in del_list:
            file_path = os.path.join(rigidsDir, f)
            if os.path.isfile(file_path):
                os.remove(file_path)
        del_list = os.listdir(elasticsDir)
        for f in del_list:
            file_path = os.path.join(elasticsDir, f)
            if os.path.isfile(file_path):
                os.remove(file_path)

        return {'FINISHED'}

class Mesh_OT_SaveSelectdMeshes2File(bpy.types.Operator):
    bl_idname = "mesh.saveselectedmeshes2file"
    bl_label = "Save selected meshes to file"
    bl_options = {'REGISTER'}

    def execute(self, context):

        basedir = os.path.dirname(bpy.data.filepath)
        if not basedir:
            raise Exception("Blend file is not saved")
        basedir += '/SPHinXsys_Blender_Addons/Mesh/'

        props.clearAllBlocksName()
        
        view_layer = bpy.context.view_layer
        obj_active = view_layer.objects.active
        allobjects = bpy.data.objects
        bpy.ops.object.select_all(action='DESELECT')

        for obj in allobjects:

            if(obj.type != 'MESH'):
                continue

            if(obj.data.m_SPHinXsysSettings.m_Export == False):
                continue

            obj.select_set(True)

            view_layer.objects.active = obj

            if(obj.data.m_SPHinXsysSettings.m_TypeofMedia == 'RigidSolid'):
                currentdir = basedir + 'Rigids/'
                props.appendBlocksName('RigidSolid', obj.name)
            elif(obj.data.m_SPHinXsysSettings.m_TypeofMedia == 'ElasticSolid'):
                currentdir = basedir + 'Elastics/'
                props.appendBlocksName('ElasticSolid', obj.name)
            elif(obj.data.m_SPHinXsysSettings.m_TypeofMedia == 'Fluid'):
                currentdir = basedir + 'Fluids/'
                props.appendBlocksName('Fluid', obj.name)

            name = bpy.path.clean_name(obj.name)
            fn = os.path.join(currentdir, name)

            bpy.ops.wm.stl_export(filepath=fn + ".stl")
            # bpy.ops.export_mesh.stl(filepath=fn + ".stl", ascii=True, use_selection=True)

            obj.select_set(False)

            print("written:", fn)
        
        # props.printBlocksName()

        view_layer.objects.active = obj_active

        return {'FINISHED'}
    
class Mesh_OT_SaveAllMeshes2File(bpy.types.Operator):
    bl_idname = "mesh.saveallmeshes2file"
    bl_label = "Save all meshes to file"
    bl_options = {'REGISTER'}

    def execute(self, context):

        basedir = os.path.dirname(bpy.data.filepath)
        if not basedir:
            raise Exception("Blend file is not saved")
        basedir += '/SPHinXsys_Blender_Addons/Mesh/'

        props.clearAllBlocksName()
        
        view_layer = bpy.context.view_layer
        obj_active = view_layer.objects.active
        allobjects = bpy.data.objects
        bpy.ops.object.select_all(action='DESELECT')

        for obj in allobjects:

            if(obj.type != 'MESH'):
                continue

            obj.select_set(True)

            view_layer.objects.active = obj

            if(obj.data.m_SPHinXsysSettings.m_TypeofMedia == 'RigidSolid'):
                currentdir = basedir + 'Rigids/'
                props.appendBlocksName('RigidSolid', obj.name)
            elif(obj.data.m_SPHinXsysSettings.m_TypeofMedia == 'ElasticSolid'):
                currentdir = basedir + 'Elastics/'
                props.appendBlocksName('ElasticSolid', obj.name)
            elif(obj.data.m_SPHinXsysSettings.m_TypeofMedia == 'Fluid'):
                currentdir = basedir + 'Fluids/'
                props.appendBlocksName('Fluid', obj.name)

            name = bpy.path.clean_name(obj.name)
            fn = os.path.join(currentdir, name)

            bpy.ops.wm.stl_export(filepath=fn + ".stl")
            #bpy.ops.export_mesh.stl(filepath=fn + ".stl", ascii=True, use_selection=True)

            obj.select_set(False)

            print("written:", fn)

        # props.printBlocksName()
        
        view_layer.objects.active = obj_active

        return {'FINISHED'}