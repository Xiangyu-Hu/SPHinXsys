import bpy
import os
import sys

path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(path)

print(sys.path)

print(sys.version)

import props

try:
    import test_3d_blender as SPHinXsys
except ImportError:
    print("SPHinXsys module failed to load!")
    isSPHinXsysInstalled = False
else:
    isSPHinXsysInstalled = True

try:
    import test_3d_blender as SPHinXsys
except ImportError as e:
    print(f"SPHinXsys module failed to load! Error: {e}")
    isSPHinXsysInstalled = False

class SPHinXsys_OT_runCase(bpy.types.Operator):
    bl_idname = "sphinxsys.runcase"
    bl_label = "Run Case"
    bl_options = {'REGISTER'}

    def execute(self, context):

        if(isSPHinXsysInstalled):
            BasicParameters = SPHinXsys.BasicParameters()
            BasicParametersFromBlender = context.scene.m_SPHinXsysGlobalSettings
            BasicParameters.sim_domain_lower = BasicParametersFromBlender.m_SimDomainLower
            BasicParameters.sim_domain_upper = BasicParametersFromBlender.m_SimDomainUpper
            BasicParameters.particle_spacing_ref = BasicParametersFromBlender.m_ParticleReferenceSpacing
            BasicParameters.gravity_g = BasicParametersFromBlender.m_Gravity
            BasicParameters.end_time = BasicParametersFromBlender.m_EndTime
            BasicParameters.rho0_f = BasicParametersFromBlender.m_Rho0f
            BasicParameters.U_ref = BasicParametersFromBlender.m_URef
            # BasicParameters.c_f = BasicParametersFromBlender.m_Cf
            BasicParameters.water_block_file_path = BasicParametersFromBlender.m_WaterBlockFilePath
            BasicParameters.rigid_block_file_path = BasicParametersFromBlender.m_RigidBlockFilePath
            BasicParameters.elastic_block_file_path = BasicParametersFromBlender.m_ElasticBlockFilePath

            w = list()
            w.append("water.stl")
            r = list()
            r.append("rigid.stl")
            e = list()
            e.append("elastic.stl")
            BasicParameters.water_block_file_name = props.m_FluidBlocksName
            BasicParameters.rigid_block_file_name = props.m_RigidBlocksName
            BasicParameters.elastic_block_file_name = props.m_ElasticBlocksName

            print(BasicParameters.water_block_file_name)
            print(BasicParameters.rigid_block_file_name)            
            print(BasicParameters.elastic_block_file_name)

            Simulator = SPHinXsys.Simulator(BasicParameters)
            Simulator.run()

            # print(BasicParameters.sim_domain_lower)
            # print(BasicParameters.sim_domain_upper)
            # print(BasicParameters.particle_spacing_ref)
            # print(BasicParameters.gravity_g)
            # print(BasicParameters.end_time)
            # print(BasicParameters.rho0_f)
            # print(BasicParameters.U_ref)
            # print(BasicParameters.c_f)
            # print(BasicParameters.water_block_file_path)
            # print(BasicParameters.rigid_block_file_path)            
            # print(BasicParameters.elastic_block_file_path)
            
        else:
            self.report({"ERROR"}, "SPHinXsys module failed to load! Run failed!")

        return {'FINISHED'}
    