import bpy

class SPHinXsys_PT_Panel(bpy.types.Panel):
    bl_idname = "SPHinXsys_PT_Panel"
    bl_label = "SPHinXsys"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "SPHinXsys"

    def draw(self, context):
        layout = self.layout

        label = layout.label(text = 'This is a label')
        #box = layout.box()
        #box.operator('object.cursor_array', text = "Save Object")
        #split = layout.split(factor=0.8)
        #split.operator('object.cursor_array', text = "Save Object")
        col = layout.column()
        col.operator('mesh.clearstlfiles', text = "Clear STL Files")
        col.operator('mesh.saveselectedmeshes2file', text = "Save Selected Meshes")
        col.operator('mesh.saveallmeshes2file', text = "Save All Meshes")
        col.operator('sphinxsys.runcase', text = 'Run Case')

        col = layout.column()
        # for prop in dir(context.scene.m_SPHinXsysGlobalSettings):
        #     if prop.startswith('m_'):
        #         col.prop(context.scene.m_SPHinXsysGlobalSettings, prop)
        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_WaterBlockFilePath')
        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_RigidBlockFilePath')
        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_ElasticBlockFilePath')

        split = layout.split(factor=0.5)
        split.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_SimDomainLower')
        split.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_SimDomainUpper')

        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_Gravity')
        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_EndTime')

        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_ParticleReferenceSpacing')
        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_Rho0f')
        col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_URef')
        # col.prop(context.scene.m_SPHinXsysGlobalSettings, 'm_Cf')
        

        meshes = bpy.data.meshes
        for item in meshes:
            box = layout.box()
            box.scale_x = 1
            box.scale_y = 1
            box.use_property_split = True
            box.use_property_decorate = False
            box.label(text = item.name)

            box.prop(item.m_SPHinXsysSettings, 'm_TypeofMedia', text = 'Media Type')
            box.prop(item.m_SPHinXsysSettings, 'm_Export', text = 'Exported?')

            if(item.m_SPHinXsysSettings.m_TypeofMedia == 'Fluid'):
                for prop in dir(item.m_SPHinXsysSettings):
                    if prop.startswith('m_Fluid'):
                        box.prop(item.m_SPHinXsysSettings, prop)
                # box.prop(item.m_SPHinXsysSettings, 'm_FluidDensity', text = 'Fluid Density')
                # box.prop(item.m_SPHinXsysSettings, 'm_FluidCharacteristicVelocity', text = 'Fluid Characteristic Velocity')
                # box.prop(item.m_SPHinXsysSettings, 'm_FluidSpeedofSound', text = 'Fluid Speed of Sound')
                # box.prop(item.m_SPHinXsysSettings, 'm_FluidReynoldsNumber', text = 'Fluid Reynolds Number')
                # box.prop(item.m_SPHinXsysSettings, 'm_FluidDynamicsViscosity', text = 'Fluid Dynamics Viscosity')
            else:
                for prop in dir(item.m_SPHinXsysSettings):
                    if prop.startswith('m_Solid'):
                        box.prop(item.m_SPHinXsysSettings, prop)
                # box.prop(item.m_SPHinXsysSettings, 'm_SolidReferenceDensity', text = 'Solid Reference Density')
                # box.prop(item.m_SPHinXsysSettings, 'm_SolidYoungsModulus', text = 'Solid Youngs Modulus')
                # box.prop(item.m_SPHinXsysSettings, 'm_SolidPoissonRatio', text = 'Solid Poisson Ratio')

            # col.label(text = item.name)
            # col.prop(item.m_SPHinXsysSettings, 'm_TypeofMedia', text = 'Media Type')
            # col.prop(item.m_SPHinXsysSettings, 'm_Export', text = 'Exported?')

            # if(item.m_SPHinXsysSettings.m_TypeofMedia == 'Fluid'):
            #     col.prop(item.m_SPHinXsysSettings, 'm_FluidDensity', text = 'Fluid Density')
            #     col.prop(item.m_SPHinXsysSettings, 'm_FluidCharacteristicVelocity', text = 'Fluid Characteristic Velocity')
            #     col.prop(item.m_SPHinXsysSettings, 'm_FluidSpeedofSound', text = 'Fluid Speed of Sound')
            #     col.prop(item.m_SPHinXsysSettings, 'm_FluidReynoldsNumber', text = 'Fluid Reynolds Number')
            #     col.prop(item.m_SPHinXsysSettings, 'm_FluidDynamicsViscosity', text = 'Fluid Dynamics Viscosity')
            # else:
            #     col.prop(item.m_SPHinXsysSettings, 'm_SolidReferenceDensity', text = 'Solid Reference Density')
            #     col.prop(item.m_SPHinXsysSettings, 'm_SolidYoungsModulus', text = 'Solid Youngs Modulus')
            #     col.prop(item.m_SPHinXsysSettings, 'm_SolidPoissonRatio', text = 'Solid Poisson Ratio')
            
